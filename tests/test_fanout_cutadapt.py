#!/usr/bin/env python3

from collections import Counter
import ctypes
import gzip
import json
import os
from pathlib import Path
import signal
import shutil
import subprocess
import struct
import sys
import tempfile
import time
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "CHARM_scripts" / "fanout_cutadapt.py"
sys.path.insert(0, str(SCRIPT.parent))
from fanout_cutadapt import parse_cutadapt_report  # noqa: E402


ADAPTER_A = "XAAAACCCCGGGG;o=12"
ADAPTER_B = "XAAAACCCCGGTT;o=12"


def audit_input_opens(directory, action):
    libc = ctypes.CDLL(None, use_errno=True)
    init = getattr(libc, "inotify_init1", None)
    add_watch = getattr(libc, "inotify_add_watch", None)
    if init is None or add_watch is None:
        raise unittest.SkipTest("Linux inotify is required")
    fd = init(os.O_NONBLOCK | os.O_CLOEXEC)
    if fd < 0:
        raise OSError(ctypes.get_errno(), "inotify_init1 failed")
    try:
        watch = add_watch(fd, os.fsencode(str(directory)), 0x00000020)
        if watch < 0:
            raise OSError(ctypes.get_errno(), "inotify_add_watch failed")
        action()
        counts = Counter()
        while True:
            try:
                data = os.read(fd, 65536)
            except BlockingIOError:
                break
            offset = 0
            while offset < len(data):
                _, mask, _, name_length = struct.unpack_from("iIII", data, offset)
                offset += struct.calcsize("iIII")
                raw_name = data[offset : offset + name_length]
                offset += name_length
                if mask & 0x00000020:
                    counts[os.fsdecode(raw_name.rstrip(b"\0"))] += 1
        return counts
    finally:
        os.close(fd)


def write_fastq(path, records, mate, bgi=False):
    with open(path, "wb") as raw_handle:
        with gzip.GzipFile(fileobj=raw_handle, mode="wb", mtime=0) as handle:
            for name, sequence in records:
                suffix = "/{} comment".format(mate) if bgi else " comment"
                record = "@{}{}\n{}\n+\n{}\n".format(
                    name, suffix, sequence, "I" * len(sequence)
                )
                handle.write(record.encode("ascii"))


def decompressed(path):
    with gzip.open(path, "rb") as handle:
        return handle.read()


def fastq_records(path):
    lines = decompressed(path).splitlines()
    if len(lines) % 4:
        raise AssertionError("truncated FASTQ: {}".format(path))
    return [tuple(lines[index : index + 4]) for index in range(0, len(lines), 4)]


def run_baseline(
    cutadapt, read1, read2, adapter, output_r1, output_r2, bgi=False, env=None
):
    raw_r1 = str(output_r1) + ".raw" if bgi else str(output_r1)
    raw_r2 = str(output_r2) + ".raw" if bgi else str(output_r2)
    command = [
        cutadapt,
        "-G",
        adapter,
        "-j",
        "1",
        "--untrimmed-output",
        os.devnull,
        "--untrimmed-paired-output",
        os.devnull,
        "-o",
        raw_r1,
        "-p",
        raw_r2,
        str(read1),
        str(read2),
    ]
    result = subprocess.run(
        command,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )
    if bgi:
        for mate, source, destination in (
            (1, raw_r1, output_r1),
            (2, raw_r2, output_r2),
        ):
            awk = subprocess.run(
                [
                    "awk",
                    '{{if(NR%4==1) gsub("/{0}"," {0}",$1); print $0}}'.format(
                        mate
                    ),
                    source,
                ],
                check=True,
                stdout=subprocess.PIPE,
            )
            with open(destination, "wb") as raw_handle:
                with gzip.GzipFile(fileobj=raw_handle, mode="wb", mtime=0) as handle:
                    handle.write(awk.stdout)
    return command, parse_cutadapt_report(result.stdout)


class FanoutCutadaptTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        environment_bin = Path(sys.executable).resolve().parent
        cls.cutadapt = shutil.which("cutadapt") or str(environment_bin / "cutadapt")
        cls.pigz = shutil.which("pigz") or str(environment_bin / "pigz")
        if not Path(cls.cutadapt).is_file():
            cls.cutadapt = None
        if not Path(cls.pigz).is_file():
            cls.pigz = None
        if not cls.cutadapt or not cls.pigz:
            raise unittest.SkipTest("cutadapt and pigz are required")

    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.read1 = self.root / "input.R1.fq.gz"
        self.read2 = self.root / "input.R2.fq.gz"
        adapter_a = ADAPTER_A.split(";", 1)[0][1:]
        adapter_b = ADAPTER_B.split(";", 1)[0][1:]
        consensus = "AAAACCCCGGGT"
        self.records_r1 = [
            ("a_only", "GATTACAGATTACA"),
            ("both", "CCCCAAAAGGGG"),
            ("b_only", "TTTTCCCCAAAA"),
            ("neither", "ACACACACACAC"),
        ]
        self.records_r2 = [
            ("a_only", adapter_a + "ACGT"),
            ("both", consensus + "TGCA"),
            ("b_only", adapter_b + "GGGG"),
            ("neither", "TGCATGCATGCAAAAA"),
        ]
        write_fastq(self.read1, self.records_r1, 1)
        write_fastq(self.read2, self.records_r2, 2)

    def tearDown(self):
        self.tempdir.cleanup()

    def helper_command(self, output_root, metrics=None, seq_format="illumina"):
        output_root.mkdir(parents=True, exist_ok=True)
        command = [
            sys.executable,
            str(SCRIPT),
            "--read1",
            str(self.read1),
            "--read2",
            str(self.read2),
            "--split",
            "a",
            ADAPTER_A,
            str(output_root / "a.R1.fq.gz"),
            str(output_root / "a.R2.fq.gz"),
            "--split",
            "b",
            ADAPTER_B,
            str(output_root / "b.R1.fq.gz"),
            str(output_root / "b.R2.fq.gz"),
            "--threads-per-split",
            "5",
            "--decompress-threads",
            "1",
            "--seq-format",
            seq_format,
            "--cutadapt",
            self.cutadapt,
            "--pigz",
            self.pigz,
            "--connect-timeout-seconds",
            "1",
            "--consumer-exit-timeout-seconds",
            "5",
        ]
        if metrics:
            command.extend(("--metrics", str(metrics)))
        return command

    def run_exact_parity(self, bgi=False):
        if bgi:
            write_fastq(self.read1, self.records_r1, 1, bgi=True)
            write_fastq(self.read2, self.records_r2, 2, bgi=True)
        baseline = self.root / "baseline"
        optimized = self.root / "optimized"
        baseline.mkdir()
        metrics_path = optimized / "metrics.json"
        baseline_counts = {}
        for name, adapter in (("a", ADAPTER_A), ("b", ADAPTER_B)):
            _, counts = run_baseline(
                self.cutadapt,
                self.read1,
                self.read2,
                adapter,
                baseline / "{}.R1.fq.gz".format(name),
                baseline / "{}.R2.fq.gz".format(name),
                bgi=bgi,
            )
            baseline_counts[name] = counts

        subprocess.run(
            self.helper_command(
                optimized, metrics=metrics_path, seq_format="bgi" if bgi else "illumina"
            ),
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        for name in ("a", "b"):
            for mate in ("R1", "R2"):
                baseline_path = baseline / "{}.{}.fq.gz".format(name, mate)
                optimized_path = optimized / "{}.{}.fq.gz".format(name, mate)
                self.assertEqual(decompressed(baseline_path), decompressed(optimized_path))
                self.assertEqual(fastq_records(baseline_path), fastq_records(optimized_path))

        metrics = json.loads(metrics_path.read_text())
        for mate in ("r1", "r2"):
            self.assertEqual(metrics["inputs"][mate]["compressed_open_count"], 1)
            self.assertEqual(
                metrics["inputs"][mate]["decompressed_bytes"],
                len(decompressed(getattr(self, "read" + mate[-1]))),
            )
        for name in ("a", "b"):
            self.assertEqual(metrics["splits"][name]["cutadapt_counts"], baseline_counts[name])
            self.assertEqual(baseline_counts[name]["input_pairs"], 4)
            self.assertEqual(baseline_counts[name]["matched_pairs"], 2)
            self.assertEqual(baseline_counts[name]["unmatched_pairs"], 2)
            self.assertEqual(baseline_counts[name]["too_short_pairs"], 0)

        observed_a = [
            record[0].split()[0][1:].decode()
            for record in fastq_records(optimized / "a.R1.fq.gz")
        ]
        observed_b = [
            record[0].split()[0][1:].decode()
            for record in fastq_records(optimized / "b.R1.fq.gz")
        ]
        self.assertEqual(observed_a, ["a_only", "both"])
        self.assertEqual(observed_b, ["both", "b_only"])

    def test_exact_parity_and_overlap_illumina(self):
        self.run_exact_parity(bgi=False)

    def test_exact_parity_bgi_header_normalization(self):
        self.run_exact_parity(bgi=True)

    def test_consumer_failure_preserves_existing_outputs(self):
        output_root = self.root / "failure"
        outputs = [
            output_root / "{}.{}.fq.gz".format(name, mate)
            for name in ("a", "b")
            for mate in ("R1", "R2")
        ]
        output_root.mkdir()
        for output in outputs:
            output.write_bytes(b"existing\n")
        command = self.helper_command(output_root)
        invalid_adapter_index = command.index(ADAPTER_B)
        command[invalid_adapter_index] = "XNOT-A-DNA-ADAPTER;o=12"
        result = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        self.assertNotEqual(result.returncode, 0)
        for output in outputs:
            self.assertEqual(output.read_bytes(), b"existing\n")
        self.assertFalse(any(output_root.glob(".multisplit.*")))

    def test_producer_failure_preserves_existing_outputs(self):
        output_root = self.root / "producer-failure"
        outputs = [
            output_root / "{}.{}.fq.gz".format(name, mate)
            for name in ("a", "b")
            for mate in ("R1", "R2")
        ]
        output_root.mkdir()
        for output in outputs:
            output.write_bytes(b"existing\n")
        self.read2.write_bytes(b"not a gzip stream\n")
        result = subprocess.run(
            self.helper_command(output_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        for output in outputs:
            self.assertEqual(output.read_bytes(), b"existing\n")
        self.assertFalse(any(output_root.glob(".multisplit.*")))

    def test_compressed_inputs_are_opened_once(self):
        baseline = self.root / "audit-baseline"
        optimized = self.root / "audit-optimized"
        baseline.mkdir()

        def run_independent_cutadapt():
            for name, adapter in (("a", ADAPTER_A), ("b", ADAPTER_B)):
                run_baseline(
                    self.cutadapt,
                    self.read1,
                    self.read2,
                    adapter,
                    baseline / "{}.R1.fq.gz".format(name),
                    baseline / "{}.R2.fq.gz".format(name),
                )

        baseline_audit = audit_input_opens(self.root, run_independent_cutadapt)

        def run_fanout():
            subprocess.run(
                self.helper_command(optimized),
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

        optimized_audit = audit_input_opens(self.root, run_fanout)
        for input_path in (self.read1, self.read2):
            self.assertEqual(baseline_audit[input_path.name], 2)
            self.assertEqual(optimized_audit[input_path.name], 1)


FAKE_PIGZ = r'''
import gzip
import os
import shutil
import sys
import time


def record_pid(role):
    directory = os.environ.get("FAKE_PID_DIR")
    if directory:
        with open(os.path.join(directory, "{}.{}".format(role, os.getpid())), "w") as handle:
            handle.write(str(os.getpid()) + "\n")


record_pid("producer")
with gzip.open(sys.argv[-1], "rb") as source:
    if os.environ.get("FAKE_PIGZ_MODE") == "stall":
        sys.stdout.buffer.write(source.read(4096))
        sys.stdout.buffer.flush()
        time.sleep(60)
    else:
        shutil.copyfileobj(source, sys.stdout.buffer, length=65536)
'''


FAKE_CUTADAPT = r'''
import gzip
import os
import shutil
import subprocess
import sys
import time


def record_pid(role, pid=None):
    pid = os.getpid() if pid is None else pid
    directory = os.environ.get("FAKE_PID_DIR")
    if directory:
        with open(os.path.join(directory, "{}.{}".format(role, pid)), "w") as handle:
            handle.write(str(pid) + "\n")


def open_output(path):
    if path.endswith(".gz"):
        return gzip.open(path, "wb")
    return open(path, "wb")


output_r1 = sys.argv[sys.argv.index("-o") + 1]
output_r2 = sys.argv[sys.argv.index("-p") + 1]
input_r1, input_r2 = sys.argv[-2:]
split_name = os.path.basename(os.path.dirname(output_r1))
mode = os.environ.get("FAKE_CUTADAPT_MODE", "normal")
record_pid("consumer")

if mode == "connect_stall":
    time.sleep(60)

with open(input_r1, "rb", buffering=0) as source_r1:
    with open(input_r2, "rb", buffering=0) as source_r2:
        if mode == "exit_one" and split_name == "b":
            sys.exit(7)
        if mode in ("never_read", "never_read_noisy", "read_then_stall"):
            child = subprocess.Popen(
                [sys.executable, "-c", "import time; time.sleep(60)"]
            )
            record_pid("descendant", child.pid)
            if mode == "read_then_stall":
                source_r1.read(4096)
                source_r2.read(4096)
            if mode == "never_read_noisy":
                while True:
                    sys.stderr.write("consumer still alive but not reading\n")
                    sys.stderr.flush()
                    time.sleep(0.01)
            time.sleep(60)
        for source, output_path in (
            (source_r1, output_r1),
            (source_r2, output_r2),
        ):
            with open_output(output_path) as destination:
                shutil.copyfileobj(source, destination, length=65536)

print("Total read pairs processed: 0")
print("Read 2 with adapter: 0")
print("Pairs written (passing filters): 0")
'''


class FanoutSupervisionTest(unittest.TestCase):
    """Deterministic process fixtures exercise backpressure and timeout cleanup."""

    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.pid_dir = self.root / "pids"
        self.pid_dir.mkdir()
        self.fake_pigz = self._write_executable("fake_pigz.py", FAKE_PIGZ)
        self.fake_cutadapt = self._write_executable(
            "fake_cutadapt.py", FAKE_CUTADAPT
        )
        self.read1 = self.root / "input.R1.fq.gz"
        self.read2 = self.root / "input.R2.fq.gz"
        write_fastq(self.read1, [("one", "ACGT")], 1)
        write_fastq(self.read2, [("one", "TGCA")], 2)

    def tearDown(self):
        self.tempdir.cleanup()

    def _write_executable(self, name, source):
        path = self.root / name
        path.write_text("#!{}\n{}".format(sys.executable, source), encoding="utf-8")
        path.chmod(0o755)
        return path

    def _outputs(self, output_root):
        return [
            output_root / "{}.{}.fq.gz".format(name, mate)
            for name in ("a", "b")
            for mate in ("R1", "R2")
        ]

    def _command(self, output_root):
        output_root.mkdir(parents=True, exist_ok=True)
        return [
            sys.executable,
            str(SCRIPT),
            "--read1",
            str(self.read1),
            "--read2",
            str(self.read2),
            "--split",
            "a",
            ADAPTER_A,
            str(output_root / "a.R1.fq.gz"),
            str(output_root / "a.R2.fq.gz"),
            "--split",
            "b",
            ADAPTER_B,
            str(output_root / "b.R1.fq.gz"),
            str(output_root / "b.R2.fq.gz"),
            "--threads-per-split",
            "1",
            "--decompress-threads",
            "1",
            "--cutadapt",
            str(self.fake_cutadapt),
            "--pigz",
            str(self.fake_pigz),
            "--connect-timeout-seconds",
            "0.35",
            "--consumer-exit-timeout-seconds",
            "0.35",
            "--stall-timeout-seconds",
            "0.35",
        ]

    def _environment(self, consumer_mode="normal", producer_mode="normal"):
        environment = os.environ.copy()
        environment.update(
            {
                "FAKE_CUTADAPT_MODE": consumer_mode,
                "FAKE_PIGZ_MODE": producer_mode,
                "FAKE_PID_DIR": str(self.pid_dir),
            }
        )
        return environment

    def _run_bounded(self, command, environment):
        started = time.monotonic()
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            env=environment,
            start_new_session=True,
        )
        try:
            stdout, stderr = process.communicate(timeout=8)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.communicate()
            self.fail("fanout test command exceeded its external 8 second guard")
        return subprocess.CompletedProcess(
            command, process.returncode, stdout, stderr
        ), time.monotonic() - started

    def _pid_is_alive(self, pid):
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return False
        return True

    def _assert_processes_reaped(self):
        pids = [int(path.read_text().strip()) for path in self.pid_dir.iterdir()]
        deadline = time.monotonic() + 2
        while time.monotonic() < deadline and any(
            self._pid_is_alive(pid) for pid in pids
        ):
            time.sleep(0.02)
        self.assertFalse(
            [pid for pid in pids if self._pid_is_alive(pid)],
            "fanout leaked fixture processes",
        )

    def _assert_failed_cleanly(
        self, consumer_mode="normal", producer_mode="normal", error_text=None
    ):
        output_root = self.root / "failure"
        outputs = self._outputs(output_root)
        output_root.mkdir()
        for output in outputs:
            output.write_bytes(b"prior-final\n")
        result, elapsed = self._run_bounded(
            self._command(output_root),
            self._environment(consumer_mode, producer_mode),
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertLess(elapsed, 4)
        if error_text:
            self.assertIn(error_text, result.stderr)
        for output in outputs:
            self.assertEqual(output.read_bytes(), b"prior-final\n")
        self.assertFalse(list(output_root.glob(".multisplit.*")))
        self._assert_processes_reaped()

    def test_consumer_alive_but_never_reads_times_out(self):
        self._large_inputs()
        self._assert_failed_cleanly(
            consumer_mode="never_read", error_text="fanout made no progress"
        )

    def test_consumer_reads_then_stalls_after_pipe_fills(self):
        self._large_inputs()
        self._assert_failed_cleanly(
            consumer_mode="read_then_stall", error_text="fanout made no progress"
        )

    def test_consumer_log_noise_does_not_mask_writer_stall(self):
        self._large_inputs()
        self._assert_failed_cleanly(
            consumer_mode="never_read_noisy",
            error_text="fanout made no progress",
        )

    def test_producer_stall_is_supervised(self):
        self._large_inputs()
        self._assert_failed_cleanly(
            producer_mode="stall", error_text="fanout made no progress"
        )

    def test_one_consumer_exit_terminates_every_child(self):
        self._large_inputs()
        self._assert_failed_cleanly(
            consumer_mode="exit_one", error_text="b consumer"
        )

    def test_fifo_connect_timeout_terminates_consumers(self):
        self._assert_failed_cleanly(
            consumer_mode="connect_stall",
            error_text="timed out waiting for Cutadapt to open FIFOs",
        )

    def test_normal_fixture_preserves_bytes_and_order(self):
        output_root = self.root / "normal"
        result, _ = self._run_bounded(
            self._command(output_root), self._environment()
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        expected = {
            "R1": decompressed(self.read1),
            "R2": decompressed(self.read2),
        }
        for output in self._outputs(output_root):
            mate = "R1" if ".R1." in output.name else "R2"
            self.assertEqual(decompressed(output), expected[mate])
        self.assertFalse(list(output_root.glob(".multisplit.*")))
        self._assert_processes_reaped()

    def _large_inputs(self):
        payload = b"ACGT" * (1024 * 1024)
        for path in (self.read1, self.read2):
            with gzip.open(path, "wb") as handle:
                handle.write(payload)


if __name__ == "__main__":
    unittest.main()

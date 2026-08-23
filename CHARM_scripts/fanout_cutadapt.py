#!/usr/bin/env python3
"""Feed one decompression of paired FASTQs to independent Cutadapt consumers."""

import argparse
import errno
import gzip
import json
import os
import re
import selectors
import signal
import shutil
import subprocess
import sys
import tempfile
import time


SPLIT_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
COUNT_PATTERNS = {
    "input_pairs": re.compile(r"Total read pairs processed:\s+([0-9,]+)"),
    "matched_pairs": re.compile(r"Read 2 with adapter:\s+([0-9,]+)"),
    "output_pairs": re.compile(r"Pairs written \(passing filters\):\s+([0-9,]+)"),
    "too_short_pairs": re.compile(r"Pairs that were too short:\s+([0-9,]+)"),
}


class SplitSpec:
    def __init__(self, name, adapter, output_r1, output_r2):
        if not SPLIT_NAME_RE.match(name):
            raise ValueError("invalid split name: {!r}".format(name))
        self.name = name
        self.adapter = adapter
        self.output_r1 = os.path.abspath(output_r1)
        self.output_r2 = os.path.abspath(output_r2)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Decompress each mate once and broadcast it to independent Cutadapt "
            "processes. A read pair may be emitted by more than one split."
        )
    )
    parser.add_argument("--read1", required=True)
    parser.add_argument("--read2", required=True)
    parser.add_argument(
        "--split",
        action="append",
        nargs=4,
        metavar=("NAME", "ADAPTER", "OUTPUT_R1", "OUTPUT_R2"),
        required=True,
    )
    parser.add_argument("--threads-per-split", type=int, required=True)
    parser.add_argument("--decompress-threads", type=int, default=1)
    parser.add_argument("--seq-format", choices=("illumina", "bgi"), default="illumina")
    parser.add_argument("--cutadapt", default="cutadapt")
    parser.add_argument("--pigz", default="pigz")
    parser.add_argument("--metrics")
    parser.add_argument("--connect-timeout-seconds", type=float, default=30.0)
    parser.add_argument("--consumer-exit-timeout-seconds", type=float, default=300.0)
    parser.add_argument(
        "--stall-timeout-seconds",
        type=float,
        help=(
            "maximum seconds without I/O or process progress; defaults to "
            "--consumer-exit-timeout-seconds for backward compatibility"
        ),
    )
    args = parser.parse_args()
    if args.threads_per_split < 1:
        parser.error("--threads-per-split must be positive")
    if args.decompress_threads < 1:
        parser.error("--decompress-threads must be positive")
    if args.connect_timeout_seconds <= 0:
        parser.error("--connect-timeout-seconds must be positive")
    if args.consumer_exit_timeout_seconds <= 0:
        parser.error("--consumer-exit-timeout-seconds must be positive")
    if args.stall_timeout_seconds is not None and args.stall_timeout_seconds <= 0:
        parser.error("--stall-timeout-seconds must be positive")
    if args.stall_timeout_seconds is None:
        args.stall_timeout_seconds = args.consumer_exit_timeout_seconds
    return args


def connect_fifos(destinations, timeout_seconds):
    deadline = time.monotonic() + timeout_seconds
    pending = [destination for values in destinations.values() for destination in values]
    while pending:
        for destination in list(pending):
            consumer = destination["consumer"]
            returncode = consumer.poll()
            if returncode is not None:
                raise RuntimeError(
                    "Cutadapt exited {} before opening {}".format(
                        returncode, destination["path"]
                    )
                )
            try:
                destination["fd"] = os.open(
                    destination["path"], os.O_WRONLY | os.O_NONBLOCK
                )
            except OSError as error:
                if error.errno != errno.ENXIO:
                    raise
            else:
                pending.remove(destination)
        if pending:
            if time.monotonic() >= deadline:
                paths = ", ".join(destination["path"] for destination in pending)
                raise RuntimeError(
                    "timed out waiting for Cutadapt to open FIFOs: {}".format(paths)
                )
            time.sleep(0.01)


def close_fifo(destination, selector=None):
    fd = destination.get("fd")
    if fd is None:
        return
    if selector is not None:
        try:
            selector.unregister(fd)
        except (KeyError, ValueError):
            pass
    try:
        os.close(fd)
    except OSError:
        pass
    destination["fd"] = None


def close_fifos(destinations, selector=None):
    for values in destinations.values():
        for destination in values:
            close_fifo(destination, selector)


def close_stream(stream, selector=None):
    if stream is None:
        return
    if selector is not None:
        try:
            selector.unregister(stream)
        except (KeyError, ValueError):
            pass
    try:
        stream.close()
    except OSError:
        pass


def signal_process_group(process, process_signal):
    try:
        os.killpg(process.pid, process_signal)
    except ProcessLookupError:
        pass
    except PermissionError:
        if process.poll() is None:
            process.send_signal(process_signal)


def terminate(processes, grace_seconds=2.0):
    """Bound cleanup and include multiprocessing descendants of each command."""
    for process in processes:
        signal_process_group(process, signal.SIGTERM)
    deadline = time.monotonic() + grace_seconds
    for process in processes:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            break
        if process.poll() is None:
            try:
                process.wait(timeout=remaining)
            except subprocess.TimeoutExpired:
                break
    for process in processes:
        # The parent may have exited while multiprocessing descendants remain.
        signal_process_group(process, signal.SIGKILL)
    for process in processes:
        if process.poll() is None:
            try:
                process.wait(timeout=grace_seconds)
            except subprocess.TimeoutExpired:
                raise RuntimeError("failed to reap child process {}".format(process.pid))


def register_capture(selector, stream, record, stream_name):
    os.set_blocking(stream.fileno(), False)
    selector.register(stream, selectors.EVENT_READ, ("capture", record, stream_name))


def register_producer_read(selector, producer):
    stream = producer["process"].stdout
    selector.register(stream, selectors.EVENT_READ, ("producer", producer))


def all_fifos_closed(producers):
    return all(producer["fifos_closed"] for producer in producers)


def process_state_summary(producers, consumers):
    parts = []
    for producer in producers:
        parts.append(
            "R{}-producer={} buffered={} eof={}".format(
                producer["mate"],
                producer["process"].poll(),
                len(producer["buffer"]),
                producer["stdout_eof"],
            )
        )
    for consumer in consumers:
        parts.append(
            "{}-consumer={}".format(
                consumer["spec"].name, consumer["process"].poll()
            )
        )
    return "; ".join(parts)


def supervise_fanout(producers, consumers, destinations, idle_timeout_seconds):
    """Move bytes through nonblocking FIFOs and supervise every child process."""
    selector = selectors.DefaultSelector()
    records = producers + consumers
    last_progress = time.monotonic()
    try:
        for producer in producers:
            stdout = producer["process"].stdout
            os.set_blocking(stdout.fileno(), False)
            register_producer_read(selector, producer)
            register_capture(selector, producer["process"].stderr, producer, "stderr")
        for consumer in consumers:
            register_capture(selector, consumer["process"].stdout, consumer, "stdout")
            register_capture(selector, consumer["process"].stderr, consumer, "stderr")

        while True:
            now = time.monotonic()
            progress = False

            for record in records:
                returncode = record["process"].poll()
                if returncode is None or record["returncode"] is not None:
                    continue
                record["returncode"] = returncode
                progress = True
                if returncode != 0:
                    raise RuntimeError(
                        "{} exited {}".format(record["label"], returncode)
                    )
                if record["kind"] == "consumer" and not all_fifos_closed(producers):
                    raise RuntimeError(
                        "{} exited before fanout completed".format(record["label"])
                    )

            if all(
                record["returncode"] is not None for record in records
            ) and not selector.get_map():
                return

            if progress:
                last_progress = time.monotonic()
            now = time.monotonic()
            remaining = idle_timeout_seconds - (now - last_progress)
            if remaining <= 0:
                raise RuntimeError(
                    "fanout made no progress for {:.3f} seconds: {}".format(
                        idle_timeout_seconds,
                        process_state_summary(producers, consumers),
                    )
                )

            events = selector.select(timeout=min(0.1, remaining))
            for key, _ in events:
                event = key.data
                if event[0] == "capture":
                    _, record, stream_name = event
                    try:
                        chunk = os.read(key.fd, 65536)
                    except BlockingIOError:
                        continue
                    if chunk:
                        record[stream_name].extend(chunk)
                    else:
                        close_stream(key.fileobj, selector)
                    # Drain diagnostics to prevent pipe backpressure, but do not
                    # let log noise conceal a stalled producer or FIFO writer.
                    continue

                if event[0] == "producer":
                    producer = event[1]
                    try:
                        chunk = os.read(key.fd, 1024 * 1024)
                    except BlockingIOError:
                        continue
                    selector.unregister(key.fileobj)
                    if not chunk:
                        close_stream(key.fileobj)
                        producer["stdout_eof"] = True
                        for destination in producer["destinations"]:
                            close_fifo(destination, selector)
                        producer["fifos_closed"] = True
                    else:
                        producer["buffer"] = chunk
                        producer["metric"]["decompressed_bytes"] += len(chunk)
                        for destination in producer["destinations"]:
                            destination["offset"] = 0
                            selector.register(
                                destination["fd"],
                                selectors.EVENT_WRITE,
                                ("destination", producer, destination),
                            )
                    progress = True
                    continue

                _, producer, destination = event
                try:
                    written = os.write(
                        destination["fd"],
                        memoryview(producer["buffer"])[destination["offset"] :],
                    )
                except BlockingIOError:
                    continue
                except BrokenPipeError:
                    raise RuntimeError(
                        "{} closed {} before fanout completed".format(
                            destination["label"], destination["path"]
                        )
                    )
                if written <= 0:
                    continue
                destination["offset"] += written
                progress = True
                if destination["offset"] == len(producer["buffer"]):
                    selector.unregister(destination["fd"])
                if all(
                    value["offset"] == len(producer["buffer"])
                    for value in producer["destinations"]
                ):
                    producer["buffer"] = b""
                    register_producer_read(selector, producer)

            if progress:
                last_progress = time.monotonic()
    finally:
        close_fifos(destinations, selector)
        for record in records:
            close_stream(record["process"].stdout, selector)
            close_stream(record["process"].stderr, selector)
        selector.close()


def parse_cutadapt_report(report):
    counts = {}
    for name, pattern in COUNT_PATTERNS.items():
        match = pattern.search(report)
        if match:
            counts[name] = int(match.group(1).replace(",", ""))
    if "input_pairs" in counts and "matched_pairs" in counts:
        counts["unmatched_pairs"] = counts["input_pairs"] - counts["matched_pairs"]
    if "too_short_pairs" not in counts:
        counts["too_short_pairs"] = 0
    return counts


def normalize_bgi(source, destination, mate):
    replacement = b" " + str(mate).encode("ascii")
    line_count = 0
    with open(source, "rb") as input_handle, open(destination, "wb") as raw_output:
        with gzip.GzipFile(fileobj=raw_output, mode="wb", mtime=0) as output_handle:
            for line_number, line in enumerate(input_handle):
                line_count = line_number + 1
                if line_number % 4 == 0:
                    fields = line.rstrip(b"\r\n").split()
                    if not fields or not fields[0].startswith(b"@"):
                        raise ValueError(
                            "malformed BGI FASTQ header at record {} in {}".format(
                                line_number // 4 + 1, source
                            )
                        )
                    fields[0] = fields[0].replace(
                        b"/" + str(mate).encode("ascii"), replacement
                    )
                    output_handle.write(b" ".join(fields) + b"\n")
                else:
                    output_handle.write(line)
    if line_count % 4 != 0:
        raise ValueError("truncated BGI FASTQ: {}".format(source))


def publish_group(staged_outputs):
    backup_root = os.path.join(os.path.dirname(staged_outputs[0][0]), "backups")
    os.makedirs(backup_root)
    backups = []
    published = []
    try:
        for index, (_, final_path) in enumerate(staged_outputs):
            os.makedirs(os.path.dirname(final_path), exist_ok=True)
            if os.path.lexists(final_path):
                backup = os.path.join(backup_root, "{:04d}".format(index))
                os.replace(final_path, backup)
                backups.append((backup, final_path))
        for staged_path, final_path in staged_outputs:
            os.replace(staged_path, final_path)
            published.append(final_path)
    except BaseException:
        for final_path in reversed(published):
            try:
                os.unlink(final_path)
            except FileNotFoundError:
                pass
        for backup, final_path in reversed(backups):
            os.replace(backup, final_path)
        raise


def main():
    args = parse_args()
    specs = [SplitSpec(*values) for values in args.split]
    if len({spec.name for spec in specs}) != len(specs):
        raise ValueError("split names must be unique")
    final_paths = [
        path for spec in specs for path in (spec.output_r1, spec.output_r2)
    ]
    if args.metrics:
        final_paths.append(os.path.abspath(args.metrics))
    if len(final_paths) != len(set(final_paths)):
        raise ValueError("all output paths must be unique")

    common_parent = os.path.commonpath([os.path.dirname(path) for path in final_paths])
    os.makedirs(common_parent, exist_ok=True)
    stage_root = tempfile.mkdtemp(prefix=".multisplit.", dir=common_parent)
    started = time.monotonic()
    all_processes = []
    destinations = {1: [], 2: []}
    consumer_records = []
    producer_records = []
    staged_outputs = []

    try:
        for spec in specs:
            split_dir = os.path.join(stage_root, spec.name)
            os.makedirs(split_dir)
            if args.seq_format == "illumina":
                cutadapt_r1 = os.path.join(split_dir, "R1.fq.gz")
                cutadapt_r2 = os.path.join(split_dir, "R2.fq.gz")
                ready_r1 = cutadapt_r1
                ready_r2 = cutadapt_r2
            else:
                cutadapt_r1 = os.path.join(split_dir, "R1.fastq")
                cutadapt_r2 = os.path.join(split_dir, "R2.fastq")
                ready_r1 = os.path.join(split_dir, "ready.R1.fq.gz")
                ready_r2 = os.path.join(split_dir, "ready.R2.fq.gz")

            # Cutadapt 2.10 multiprocessing can reopen named FIFOs reliably;
            # inherited /dev/fd paths hang when more than one core is used.
            input_r1 = os.path.join(split_dir, "input.R1.fastq")
            input_r2 = os.path.join(split_dir, "input.R2.fastq")
            os.mkfifo(input_r1)
            os.mkfifo(input_r2)
            command = [
                args.cutadapt,
                "-G",
                spec.adapter,
                "-j",
                str(args.threads_per_split),
                "--untrimmed-output",
                os.devnull,
                "--untrimmed-paired-output",
                os.devnull,
                "-o",
                cutadapt_r1,
                "-p",
                cutadapt_r2,
                input_r1,
                input_r2,
            ]
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                start_new_session=True,
            )
            all_processes.append(process)
            consumer = {
                "kind": "consumer",
                "label": "{} consumer".format(spec.name),
                "spec": spec,
                "process": process,
                "returncode": None,
                "stdout": bytearray(),
                "stderr": bytearray(),
            }
            consumer_records.append(consumer)
            destinations[1].append(
                {
                    "path": input_r1,
                    "consumer": process,
                    "label": consumer["label"],
                    "fd": None,
                    "offset": 0,
                }
            )
            destinations[2].append(
                {
                    "path": input_r2,
                    "consumer": process,
                    "label": consumer["label"],
                    "fd": None,
                    "offset": 0,
                }
            )
            staged_outputs.extend(
                ((ready_r1, spec.output_r1), (ready_r2, spec.output_r2))
            )

        # Producers start only after every FIFO reader is connected. All FIFO
        # descriptors remain nonblocking for the lifetime of the fanout.
        connect_fifos(destinations, args.connect_timeout_seconds)

        input_metrics = {}
        for mate, input_path in ((1, args.read1), (2, args.read2)):
            command = [
                args.pigz,
                "-dc",
                "-p",
                str(args.decompress_threads),
                input_path,
            ]
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                start_new_session=True,
            )
            all_processes.append(process)
            metric = {
                "path": os.path.abspath(input_path),
                "compressed_bytes": os.path.getsize(input_path),
                "compressed_open_count": 1,
                "decompressed_bytes": 0,
            }
            input_metrics["r{}".format(mate)] = metric
            producer_records.append(
                {
                    "kind": "producer",
                    "label": "R{} producer".format(mate),
                    "mate": mate,
                    "process": process,
                    "returncode": None,
                    "stderr": bytearray(),
                    "metric": metric,
                    "destinations": destinations[mate],
                    "buffer": b"",
                    "stdout_eof": False,
                    "fifos_closed": False,
                }
            )

        supervise_fanout(
            producer_records,
            consumer_records,
            destinations,
            args.stall_timeout_seconds,
        )

        split_metrics = {}
        for consumer in consumer_records:
            spec = consumer["spec"]
            stdout = bytes(consumer["stdout"]).decode("utf-8", "replace")
            stderr = bytes(consumer["stderr"]).decode("utf-8", "replace")
            if stdout:
                sys.stderr.write("[{} cutadapt stdout]\n{}".format(spec.name, stdout))
                if not stdout.endswith("\n"):
                    sys.stderr.write("\n")
            if stderr:
                sys.stderr.write("[{} cutadapt stderr]\n{}".format(spec.name, stderr))
                if not stderr.endswith("\n"):
                    sys.stderr.write("\n")
            split_metrics[spec.name] = {
                "adapter": spec.adapter,
                "cutadapt_counts": parse_cutadapt_report(stdout),
                "returncode": consumer["returncode"],
            }

        if args.seq_format == "bgi":
            for spec in specs:
                split_dir = os.path.join(stage_root, spec.name)
                normalize_bgi(
                    os.path.join(split_dir, "R1.fastq"),
                    os.path.join(split_dir, "ready.R1.fq.gz"),
                    1,
                )
                normalize_bgi(
                    os.path.join(split_dir, "R2.fastq"),
                    os.path.join(split_dir, "ready.R2.fq.gz"),
                    2,
                )

        metrics = {
            "contract": "independent_cutadapt_consumers_v1",
            "elapsed_seconds": time.monotonic() - started,
            "inputs": input_metrics,
            "splits": split_metrics,
        }
        if args.metrics:
            staged_metrics = os.path.join(stage_root, "metrics.json")
            with open(staged_metrics, "w", encoding="ascii") as handle:
                json.dump(metrics, handle, indent=2, sort_keys=True)
                handle.write("\n")
            staged_outputs.append((staged_metrics, os.path.abspath(args.metrics)))

        publish_group(staged_outputs)
        sys.stderr.write(
            "fanout metrics: {}\n".format(json.dumps(metrics, sort_keys=True))
        )
    except BaseException:
        close_fifos(destinations)
        terminate(all_processes)
        raise
    finally:
        for process in all_processes:
            close_stream(process.stdout)
            close_stream(process.stderr)
        shutil.rmtree(stage_root, ignore_errors=True)


if __name__ == "__main__":
    main()

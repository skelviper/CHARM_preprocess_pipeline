#!/usr/bin/env python3

import csv
import gzip
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam


PIPELINE = Path(__file__).resolve().parents[1]
SCRIPTS = PIPELINE / "CHARM_scripts"
sys.path.insert(0, str(SCRIPTS))

from r2_fragment_selector_core import cigar_metrics, select_fragment  # noqa: E402


HEADER_DICT = {
    "HD": {"VN": "1.6"},
    "SQ": [{"SN": "chr1", "LN": 10000}],
}


def make_record(qname, start, cigar, origin="R2", sequence=None, extra_flag=0, mapq=60):
    header = pysam.AlignmentHeader.from_dict(HEADER_DICT)
    record = pysam.AlignedSegment(header)
    record.query_name = qname
    flag = 1 | (64 if origin == "R1" else 128 if origin == "R2" else 0)
    record.flag = flag | extra_flag
    record.reference_id = 0
    record.reference_start = start
    record.mapping_quality = mapq
    record.cigartuples = cigar
    sequence_length = sum(
        length for operation, length in cigar if operation in (0, 1, 4)
    )
    record.query_sequence = sequence or ("T" * sequence_length)
    record.query_qualities = pysam.qualitystring_to_array("I" * sequence_length)
    return record


def read_gzip_bytes(path):
    with gzip.open(str(path), "rb") as handle:
        return handle.read()


class ContractUnitTests(unittest.TestCase):
    def test_terminal_clip_sums_soft_and_hard_clips(self):
        metrics = cigar_metrics([(5, 2), (4, 3), (0, 10), (4, 4), (5, 5)])
        self.assertEqual(metrics["left_clip"], 5)
        self.assertEqual(metrics["right_clip"], 9)
        self.assertEqual(metrics["aligned_query_length"], 10)
        self.assertEqual(metrics["full_query_length"], 24)

    def test_selection_uses_qstart_then_primary(self):
        base = {
            "qname": "q",
            "origin": "R2",
            "chrom": "chr1",
            "end": 110,
            "strand": "-",
            "mapq": 60,
            "q_end": 10,
            "aligned_query_length": 10,
            "phase": ".",
            "phase0_count": 0,
            "phase1_count": 0,
            "other_allele_count": 0,
        }
        clipped_primary = dict(base, flag=129, start=100, q_start=2)
        unclipped_supplementary = dict(base, flag=129 | 2048, start=200, q_start=0)
        result = select_fragment([clipped_primary, unclipped_supplementary])
        self.assertEqual(result["selected"]["start"], 200)

        primary = dict(base, flag=129, start=300, q_start=0)
        result = select_fragment([unclipped_supplementary, primary])
        self.assertEqual(result["selected"]["start"], 300)

    def test_direct_r2_phase_precedes_r1_fallback(self):
        base = {
            "qname": "q",
            "chrom": "chr1",
            "end": 110,
            "strand": "-",
            "mapq": 60,
            "q_start": 0,
            "q_end": 10,
            "aligned_query_length": 10,
            "phase0_count": 1,
            "phase1_count": 0,
            "other_allele_count": 0,
        }
        r2 = dict(base, origin="R2", flag=129, start=100, phase="0")
        r1 = dict(base, origin="R1", flag=65, start=200, phase="1")
        result = select_fragment([r1, r2])
        self.assertEqual(result["final_phase"], "0")
        self.assertEqual(result["phase_source"], "R2")
        self.assertTrue(result["r2_r1_conflict"])


class BackendParityTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.work = Path(self.temporary.name)
        self.bam = self.work / "synthetic.bam"
        self.snp = self.work / "phased_snp.tsv"
        self.index_manifest = self.work / "snp_index" / "manifest.tsv"

        records = [
            make_record("singleton_r2", 100, [(0, 10)], sequence="A" + "T" * 9),
            make_record("r1_fallback", 200, [(0, 10)], origin="R1", sequence="G" + "T" * 9),
            make_record("r1_fallback", 300, [(0, 10)]),
            make_record("clipped_r2", 400, [(4, 2), (0, 8)]),
            make_record("supplementary_wins", 500, [(4, 2), (0, 8)]),
            make_record(
                "supplementary_wins",
                600,
                [(0, 5), (5, 5)],
                extra_flag=2048,
                sequence="TTTTT",
            ),
            make_record("secondary_only", 650, [(0, 10)], extra_flag=256),
            make_record(
                "reverse_accept",
                700,
                [(4, 2), (0, 8)],
                sequence="TT" + "G" + "T" * 7,
                extra_flag=16,
            ),
            make_record(
                "reverse_reject", 750, [(0, 8), (4, 2)], extra_flag=16
            ),
            make_record(
                "r2_direct_precedence",
                800,
                [(0, 10)],
                origin="R1",
                sequence="G" + "T" * 9,
            ),
            make_record(
                "r2_direct_precedence",
                900,
                [(0, 10)],
                sequence="A" + "T" * 9,
            ),
            make_record(
                "deletion_then_match",
                950,
                [(0, 3), (2, 2), (0, 5)],
                sequence="TTTGTTTT",
            ),
            make_record(
                "insertion_then_match",
                1000,
                [(0, 3), (1, 2), (0, 5)],
                sequence="TTTTTATTTT",
            ),
        ]
        with pysam.AlignmentFile(str(self.bam), "wb", header=HEADER_DICT) as handle:
            for record in reversed(records):
                handle.write(record)

        self.snp.write_text(
            "chr1\t101\tA\tG\n"
            "chr1\t201\tA\tG\n"
            "chr1\t701\tA\tG\n"
            "chr1\t801\tA\tG\n"
            "chr1\t901\tA\tG\n"
            "chr1\t956\tA\tG\n"
            "chr1\t1004\tA\tG\n"
        )
        subprocess.run(
            [
                sys.executable,
                str(SCRIPTS / "build_phased_snp_index.py"),
                "--input",
                str(self.snp),
                "--output",
                str(self.index_manifest),
            ],
            check=True,
        )

    def tearDown(self):
        self.temporary.cleanup()

    def run_backend(self, backend, label=None, batch_qnames=None):
        output_dir = self.work / (label or backend)
        output_dir.mkdir()
        command = [
            sys.executable,
            str(SCRIPTS / "select_r2_5prime_unclipped.py"),
            "--sample",
            "synthetic",
            "--split",
            "atac",
            "--bam",
            str(self.bam),
            "--backend",
            backend,
            "--output",
            str(output_dir / "fragments.bed.gz"),
            "--phase-output",
            str(output_dir / "phase.tsv.gz"),
            "--stats",
            str(output_dir / "stats.tsv"),
            "--log",
            str(output_dir / "run.log"),
            "--temp-root",
            str(output_dir / "tmp"),
            "--threads",
            "1",
        ]
        if backend == "python":
            command.extend(["--snp-index", str(self.index_manifest)])
            if batch_qnames is not None:
                command.extend(["--python-batch-qnames", str(batch_qnames)])
        else:
            command.extend(
                [
                    "--hickit",
                    str(SCRIPTS / "hickit.js"),
                    "--snp",
                    str(self.snp),
                ]
            )
        subprocess.run(command, check=True)
        return output_dir

    def test_python_and_hickit_reference_are_identical(self):
        python_output = self.run_backend("python")
        python_batch_one = self.run_backend(
            "python", label="python_batch_one", batch_qnames=1
        )
        hickit_output = self.run_backend("hickit_reference")
        for filename in ("fragments.bed.gz", "phase.tsv.gz"):
            self.assertEqual(
                read_gzip_bytes(python_output / filename),
                read_gzip_bytes(hickit_output / filename),
                filename,
            )
            self.assertEqual(
                read_gzip_bytes(python_output / filename),
                read_gzip_bytes(python_batch_one / filename),
                "{} batch-size invariance".format(filename),
            )
        self.assertEqual(
            (python_output / "stats.tsv").read_bytes(),
            (hickit_output / "stats.tsv").read_bytes(),
        )
        self.assertEqual(
            (python_output / "stats.tsv").read_bytes(),
            (python_batch_one / "stats.tsv").read_bytes(),
        )

        bed_rows = read_gzip_bytes(python_output / "fragments.bed.gz").decode().splitlines()
        self.assertEqual(len(bed_rows), 7)
        by_qname = {row.split("\t")[3]: row.split("\t") for row in bed_rows}
        self.assertEqual(by_qname["supplementary_wins"][1:3], ["600", "605"])
        self.assertEqual(by_qname["reverse_accept"][5], "+")
        self.assertNotIn("clipped_r2", by_qname)
        self.assertNotIn("secondary_only", by_qname)
        self.assertNotIn("reverse_reject", by_qname)

        phase_text = read_gzip_bytes(python_output / "phase.tsv.gz").decode()
        phase_rows = list(csv.DictReader(phase_text.splitlines(), delimiter="\t"))
        phase_by_qname = {row["qname"]: row for row in phase_rows}
        self.assertEqual(phase_by_qname["singleton_r2"]["phase_source"], "R2")
        self.assertEqual(phase_by_qname["r1_fallback"]["phase_source"], "R1")
        self.assertEqual(
            phase_by_qname["r2_direct_precedence"]["final_phase"], "0"
        )
        self.assertEqual(
            phase_by_qname["r2_direct_precedence"]["r2_r1_conflict"], "1"
        )
        self.assertEqual(phase_by_qname["deletion_then_match"]["final_phase"], "1")
        self.assertEqual(phase_by_qname["insertion_then_match"]["final_phase"], "0")


if __name__ == "__main__":
    unittest.main()

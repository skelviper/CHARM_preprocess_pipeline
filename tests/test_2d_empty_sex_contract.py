#!/usr/bin/env python3

import csv
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
CLASSIFIER = ROOT / "CHARM_scripts" / "classify_2d_bam.py"
OUTPUT_UTILS = ROOT / "CHARM_scripts" / "output_utils.sh"


def sam_record(qname, chromosome="chr1", mapq=60, flag=0):
    if flag & 0x4:
        return "{}\t{}\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n".format(
            qname, flag
        )
    return "{}\t{}\t{}\t1\t{}\t4M\t*\t0\t0\tACGT\tIIII\n".format(
        qname, flag, chromosome, mapq
    )


def read_single_tsv(path):
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise AssertionError("expected one row: {}".format(rows))
    return rows[0]


@unittest.skipUnless(shutil.which("samtools"), "samtools is not on PATH")
class BamSexContractTests(unittest.TestCase):
    def make_bam(self, root, name, records):
        sam = root / (name + ".sam")
        bam = root / (name + ".bam")
        sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:1000\n"
            "@SQ\tSN:chrX\tLN:1000\n"
            "@SQ\tSN:chrY\tLN:1000\n"
            + "".join(records)
        )
        completed = subprocess.run(
            ["samtools", "view", "-b", "-o", str(bam), str(sam)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        self.assertEqual(0, completed.returncode, completed.stderr)
        return bam

    def classify(self, root, bam, fallback="XX", threshold="0.015"):
        classification = root / (bam.stem + ".sex.tsv")
        yperx = root / (bam.stem + ".yperx.txt")
        completed = subprocess.run(
            [
                sys.executable,
                str(CLASSIFIER),
                "--bam",
                str(bam),
                "--classification-output",
                str(classification),
                "--yperx-output",
                str(yperx),
                "--yperx-threshold",
                threshold,
                "--no-xy-fallback",
                fallback,
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        self.assertEqual(0, completed.returncode, completed.stderr)
        return read_single_tsv(classification), yperx.read_text()

    def test_empty_alignment_states_are_distinct(self):
        cases = (
            ("header", [], "header_only", "0", "0"),
            ("unmapped", [sam_record("u1", flag=4)], "unmapped_only", "1", "0"),
            ("low", [sam_record("l1", mapq=19)], "low_mapq_only", "1", "1"),
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name, records, state, total, mapped in cases:
                with self.subTest(state=state):
                    bam = self.make_bam(root, name, records)
                    row, yperx = self.classify(root, bam)
                    self.assertEqual(state, row["input_state"])
                    self.assertEqual("not_applicable", row["sex_state"])
                    self.assertEqual("false", row["measured"])
                    self.assertEqual(total, row["total_records"])
                    self.assertEqual(mapped, row["mapped_records"])
                    self.assertEqual("0", row["qualifying_records"])
                    self.assertEqual("", yperx)

    def test_zero_denominator_states_and_fallback_are_explicit(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            autosomal = self.make_bam(
                root, "autosomal", [sam_record("a1"), sam_record("a2")]
            )
            xx, xx_scalar = self.classify(root, autosomal, fallback="XX")
            xy, xy_scalar = self.classify(root, autosomal, fallback="XY")
            for row in (xx, xy):
                self.assertEqual("undetermined_no_xy", row["sex_state"])
                self.assertEqual("undefined_zero_zero", row["ratio_state"])
                self.assertEqual("false", row["measured"])
                self.assertEqual("0", row["denominator_x"])
                self.assertEqual("", row["measured_yperx"])
            self.assertEqual("XX", xx["sex_call"])
            self.assertEqual("XY", xy["sex_call"])
            self.assertEqual("0.000000\n", xx_scalar)
            self.assertNotEqual(xx_scalar, xy_scalar)

            y_only = self.make_bam(root, "y_only", [sam_record("y1", "chrY")])
            row, scalar = self.classify(root, y_only)
            self.assertEqual("xy_y_only", row["sex_state"])
            self.assertEqual("XY", row["sex_call"])
            self.assertEqual("positive_infinity", row["ratio_state"])
            self.assertEqual("false", row["measured"])
            self.assertEqual("0", row["denominator_x"])
            self.assertEqual("1.000000\n", scalar)

            x_only = self.make_bam(root, "x_only", [sam_record("x1", "chrX")])
            row, scalar = self.classify(root, x_only)
            self.assertEqual("xx_x_only", row["sex_state"])
            self.assertEqual("XX", row["sex_call"])
            self.assertEqual("true", row["measured"])
            self.assertEqual("0.000000\n", scalar)

    def test_ratio_threshold_is_strict_and_uses_six_decimal_legacy_value(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name, y_count, expected_state, expected_call, expected_ratio in (
                ("below", 14, "xx_ratio", "XX", "0.014000\n"),
                ("above", 16, "xy_ratio", "XY", "0.016000\n"),
            ):
                records = [sam_record("x{}".format(i), "chrX") for i in range(1000)]
                records += [sam_record("y{}".format(i), "chrY") for i in range(y_count)]
                bam = self.make_bam(root, name, records)
                row, scalar = self.classify(root, bam)
                self.assertEqual(expected_state, row["sex_state"])
                self.assertEqual(expected_call, row["sex_call"])
                self.assertEqual("true", row["measured"])
                self.assertEqual(expected_ratio, scalar)

    def test_empty_artifacts_are_valid_and_classifier_failures_are_nonzero(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name, records, expected_state in (
                ("header", [], "header_only"),
                ("unmapped", [sam_record("u1", flag=4)], "unmapped_only"),
                ("low", [sam_record("l1", mapq=19)], "low_mapq_only"),
            ):
                with self.subTest(state=expected_state):
                    bam = self.make_bam(root, name, records)
                    final = root / (name + ".final")
                    environment = os.environ.copy()
                    environment.update(
                        {
                            "BAM": str(bam),
                            "CLASSIFIER": str(CLASSIFIER),
                            "OUTPUT_UTILS": str(OUTPUT_UTILS),
                            "FINAL": str(final),
                            "STATE": expected_state,
                        }
                    )
                    completed = subprocess.run(
                        ["bash", "-c", r'''
set -euo pipefail
source "$OUTPUT_UTILS"
mkdir -p "$FINAL"
summary=$(python "$CLASSIFIER" --bam "$BAM" \
  --classification-output "$FINAL/sex.tsv" \
  --yperx-output "$FINAL/yperx.txt" \
  --yperx-threshold 0.015)
IFS='|' read -r input_state sex_state sex_call qualifying <<< "$summary"
test "$input_state" = "$STATE"
test "$qualifying" -eq 0
charm_make_empty_seg_from_bam "$BAM" "$FINAL/contacts.seg.gz"
test "$(charm_gzip_nonheader_records "$FINAL/contacts.seg.gz")" -eq 0
'''],
                        env=environment,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        check=False,
                    )
                    self.assertEqual(0, completed.returncode, completed.stderr)
                    self.assertEqual("", (final / "yperx.txt").read_text())
                    self.assertEqual(
                        expected_state,
                        read_single_tsv(final / "sex.tsv")["input_state"],
                    )
                    with gzip.open(final / "contacts.seg.gz", "rt") as handle:
                        self.assertTrue(
                            all(line.startswith("#chromosome:") for line in handle)
                        )

            fake = root / "fake-samtools"
            fake.write_text(
                "#!/usr/bin/env bash\n"
                "printf 'partial\\t0\\tchr1\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\tIIII\\n'\n"
                "printf 'injected failure\\n' >&2\n"
                "exit 42\n"
            )
            fake.chmod(0o755)
            failed = root / "failed"
            failed.mkdir()
            environment = os.environ.copy()
            environment.update(
                {
                    "CLASSIFIER": str(CLASSIFIER),
                    "FAKE": str(fake),
                    "FINAL": str(failed),
                }
            )
            completed = subprocess.run(
                ["bash", "-c", r'''
set -euo pipefail
python "$CLASSIFIER" --bam ignored.bam --samtools "$FAKE" \
  --classification-output "$FINAL/sex.tsv" \
  --yperx-output "$FINAL/yperx.txt" \
  --yperx-threshold 0.015
'''],
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            self.assertNotEqual(0, completed.returncode)


if __name__ == "__main__":
    unittest.main()

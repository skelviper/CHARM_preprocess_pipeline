#!/usr/bin/env python3
"""Focused tests for record counting and valid empty-output contracts."""

import gzip
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_UTILS = REPO_ROOT / "CHARM_scripts" / "output_utils.sh"


class OutputUtilityTests(unittest.TestCase):
    def run_bash(self, script, environment):
        env = os.environ.copy()
        env.update({key: str(value) for key, value in environment.items()})
        return subprocess.run(
            ["bash", "-c", script],
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

    def test_biological_empty_pairs_are_format_valid(self):
        with tempfile.TemporaryDirectory(prefix="charm-output-empty-pairs-") as temporary:
            root = Path(temporary)
            seg = root / "empty.seg.gz"
            with gzip.open(seg, "wt") as handle:
                handle.write("#chromosome: chr1 1000\n")
                handle.write("#chromosome: chrX 500\n")
            output = root / "empty.pairs.gz"
            completed = self.run_bash(
                r"""
set -euo pipefail
source "$OUTPUT_UTILS"
charm_make_empty_pairs_from_seg "$SEG" "$OUTPUT" gzip
test "$(charm_gzip_nonheader_records "$OUTPUT")" -eq 0
""",
                {"OUTPUT_UTILS": OUTPUT_UTILS, "SEG": seg, "OUTPUT": output},
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with gzip.open(output, "rt") as handle:
                lines = handle.read().splitlines()
            self.assertEqual(lines[0], "## pairs format v1.0")
            self.assertIn("#chromosome: chr1 1000", lines)
            self.assertEqual(
                lines[-1],
                "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1",
            )

    @unittest.skipUnless(shutil.which("samtools"), "samtools is not on PATH")
    def test_header_only_bam_produces_format_valid_empty_seg(self):
        with tempfile.TemporaryDirectory(prefix="charm-output-empty-bam-") as temporary:
            root = Path(temporary)
            bam = root / "empty.bam"
            seg = root / "empty.seg.gz"
            completed = self.run_bash(
                r"""
set -euo pipefail
source "$OUTPUT_UTILS"
printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chrX\tLN:500\n' | \
  samtools view -b -o "$BAM" -
test "$(samtools view -c "$BAM")" -eq 0
charm_make_empty_seg_from_bam "$BAM" "$SEG"
test "$(charm_gzip_nonheader_records "$SEG")" -eq 0
""",
                {"OUTPUT_UTILS": OUTPUT_UTILS, "BAM": bam, "SEG": seg},
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with gzip.open(seg, "rt") as handle:
                self.assertEqual(
                    handle.read().splitlines(),
                    ["#chromosome: chr1 1000", "#chromosome: chrX 500"],
                )

    @unittest.skipUnless(shutil.which("bgzip"), "bgzip is not on PATH")
    def test_biological_empty_impute_pairs_are_valid_bgzip(self):
        with tempfile.TemporaryDirectory(prefix="charm-output-empty-bgzip-") as temporary:
            root = Path(temporary)
            source = root / "source.pairs.gz"
            with gzip.open(source, "wt") as handle:
                handle.write("## pairs format v1.0\n")
                handle.write("#sorted: chr1-chr2-pos1-pos2\n")
                handle.write("#shape: upper triangle\n")
                handle.write("#chromosome: chr1 1000\n")
                handle.write(
                    "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1\n"
                )
            output = root / "empty.impute.pairs.gz"
            completed = self.run_bash(
                r"""
set -euo pipefail
source "$OUTPUT_UTILS"
charm_copy_gzip_headers "$SOURCE" "$OUTPUT" bgzip
test "$(charm_gzip_nonheader_records "$OUTPUT")" -eq 0
""",
                {
                    "OUTPUT_UTILS": OUTPUT_UTILS,
                    "SOURCE": source,
                    "OUTPUT": output,
                },
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with gzip.open(source, "rt") as source_handle:
                expected = source_handle.read()
            with gzip.open(output, "rt") as output_handle:
                self.assertEqual(output_handle.read(), expected)

    @unittest.skipUnless(
        shutil.which("bgzip") and shutil.which("tabix"),
        "bgzip and tabix are not on PATH",
    )
    def test_empty_combined_fragment_bgzip_has_valid_index(self):
        with tempfile.TemporaryDirectory(prefix="charm-output-empty-tabix-") as temporary:
            output = Path(temporary) / "empty.fragments.bgz"
            completed = self.run_bash(
                r"""
set -euo pipefail
bgzip -c /dev/null > "$OUTPUT"
tabix -f -p bed "$OUTPUT"
test -s "$OUTPUT.tbi"
""",
                {"OUTPUT": output},
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)


class RuleAuditTests(unittest.TestCase):
    RULES = (
        REPO_ROOT / "rules" / "CHARM_split.rules",
        REPO_ROOT / "rules" / "CHARM_cuttag.rules",
        REPO_ROOT / "rules" / "scHiC_2dprocess.rules",
        REPO_ROOT / "rules" / "scHiC_3dprocess.rules",
    )

    def test_lifecycle_publication_framework_is_absent(self):
        text = "\n".join(path.read_text() for path in self.RULES)
        helper = OUTPUT_UTILS.read_text()
        for forbidden in (
            "charm_fc_",
            "CHARM_FC_",
            "fail_closed_io.sh",
            ".status.tsv",
            "charm_fc_begin",
            "charm_fc_finish",
            "charm_validate_",
            "quickcheck",
            "idxstats",
            "tabix -l",
            "test -f",
        ):
            self.assertNotIn(forbidden, text)
            self.assertNotIn(forbidden, helper)
        for forbidden in ("trap ", "mktemp", "mv --"):
            self.assertNotIn(forbidden, helper)

    def test_known_failure_masking_patterns_are_absent(self):
        text = "\n".join(path.read_text() for path in self.RULES)
        self.assertNotIn("writing empty", text)
        self.assertNotIn("|| touch", text)
        self.assertNotIn("if !", text)

    def test_scientific_commands_and_thresholds_are_retained(self):
        two_d = (REPO_ROOT / "rules" / "scHiC_2dprocess.rules").read_text()
        three_d = (REPO_ROOT / "rules" / "scHiC_3dprocess.rules").read_text()
        cuttag = (REPO_ROOT / "rules" / "CHARM_cuttag.rules").read_text()
        for token in ("-5SP", "--dup-dist=0", "--dup-dist=500", "--dup-dist=1"):
            self.assertIn(token, two_d)
        for token in ("clean_leg", "clean_isolated", "clean_splicing"):
            self.assertIn(token, two_d)
        for token in ("-Sr1m", "-b200k", "-D5 -b50k", "-D5 -b20k"):
            self.assertIn(token, three_d)
        self.assertIn("--min-mapq {params.min_mapq}", cuttag)
        self.assertIn("--min-baseq {params.min_baseq}", cuttag)

    def test_changed_rules_use_shell_strict_mode_and_only_needed_helpers(self):
        strict_expected = {
            "CHARM_split.rules": ["multiSplit"],
            "CHARM_cuttag.rules": [
                "R2_mapping",
                "bam2bed",
                "generate_dedup_fragments",
                "foramt_fragmenets",
                "combine_fragments",
            ],
            "scHiC_2dprocess.rules": [
                "bwa_map",
                "align2pairs",
                "seg2pairs",
                "clean_leg",
                "clean_isolated",
                "clean_splicing",
                "hickit_2d",
            ],
            "scHiC_3dprocess.rules": [
                "hickit_3d",
                "compress_3dg",
                "hickit_clean3D",
                "compress_clean3d",
                "hickit_align3D",
            ],
        }
        helper_expected = {
            "scHiC_2dprocess.rules": {
                "align2pairs",
                "seg2pairs",
                "clean_leg",
                "clean_isolated",
                "clean_splicing",
                "hickit_2d",
            },
            "scHiC_3dprocess.rules": {
                "hickit_3d",
                "hickit_clean3D",
                "hickit_align3D",
            },
        }
        for filename, rule_names in strict_expected.items():
            lines = (REPO_ROOT / "rules" / filename).read_text().splitlines()
            starts = [
                (number, line.split()[1].rstrip(":"))
                for number, line in enumerate(lines)
                if line.startswith("rule ")
            ]
            for index, (start, name) in enumerate(starts):
                if name not in rule_names:
                    continue
                end = starts[index + 1][0] if index + 1 < len(starts) else len(lines)
                block = "\n".join(lines[start:end])
                self.assertIn("set -euo pipefail", block, name)
                if name in helper_expected.get(filename, set()):
                    self.assertIn("output_utils.sh", block, name)
                    self.assertIn("source {input.output_utils}", block, name)
                else:
                    self.assertNotIn("output_utils.sh", block, name)
                    self.assertNotIn("source {input.output_utils}", block, name)


@unittest.skipUnless(shutil.which("snakemake"), "snakemake is not on PATH")
class SnakemakeOutputCleanupTests(unittest.TestCase):
    def test_failed_rule_removes_all_declared_partial_outputs(self):
        with tempfile.TemporaryDirectory(prefix="charm-snakemake-cleanup-") as temporary:
            root = Path(temporary)
            (root / "Snakefile").write_text(
                "rule all:\n"
                "    input: 'one.txt', 'two.txt'\n\n"
                "rule broken:\n"
                "    output: 'one.txt', 'two.txt'\n"
                "    shell:\n"
                "        \"printf one > {output[0]}; "
                "printf two > {output[1]}; exit 17\"\n"
            )
            completed = subprocess.run(
                ["snakemake", "--cores", "1"],
                cwd=str(root),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
            self.assertNotEqual(completed.returncode, 0, completed.stdout)
            self.assertFalse((root / "one.txt").exists(), completed.stdout)
            self.assertFalse((root / "two.txt").exists(), completed.stdout)


if __name__ == "__main__":
    unittest.main()

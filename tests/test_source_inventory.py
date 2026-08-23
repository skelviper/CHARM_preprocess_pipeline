#!/usr/bin/env python3

import csv
import os
from pathlib import Path
import sys
import tempfile
import time
import unittest


PIPELINE = Path(__file__).resolve().parents[1]
SCRIPTS = PIPELINE / "CHARM_scripts"
sys.path.insert(0, str(SCRIPTS))

from source_inventory import (  # noqa: E402
    SourceInventoryError,
    build_source_inventory,
    write_source_hash_inventory,
)


class SourceInventoryTests(unittest.TestCase):
    def test_inventory_covers_execution_surface_and_excludes_tests(self):
        entries = build_source_inventory(
            str(PIPELINE), [str(PIPELINE / "config.yaml")]
        )
        relative = {
            Path(os.path.relpath(Path(path).resolve(), PIPELINE.resolve())).as_posix()
            for _, path in entries
        }
        required = {
            "runCHARM.sh",
            "CHARM.smk",
            "rules/scHiC_3dprocess.rules",
            "CHARM_scripts/input_contract.py",
            "CHARM_scripts/classify_2d_bam.py",
            "CHARM_scripts/classify_rna_gene_compatibility.py",
            "CHARM_scripts/audit_complete_run.py",
            "CHARM_scripts/generate_stat_contract.py",
            "CHARM_scripts/resolve_run_config.py",
            "CHARM_scripts/source_inventory.py",
            "CHARM_scripts/summarize_rna_output_modes.py",
            "envs/charm.yml",
            "../CHARMtools/UPSTREAM_MANIFEST.sha256",
            "../CHARMtools/charm_preprocess/clean_leg.py",
        }
        self.assertTrue(required.issubset(relative))
        self.assertFalse(any(path.startswith("tests/") for path in relative))
        self.assertFalse(any(path.startswith("docs/") for path in relative))
        self.assertFalse(any(path.startswith("profiles/") for path in relative))
        self.assertFalse(any(path.startswith("vendor/") for path in relative))
        self.assertFalse(
            any(path.startswith("../CHARMtools/ref/M23_cpg/") for path in relative)
        )

    def test_inventory_and_hash_output_are_deterministic_and_content_stable(self):
        first = build_source_inventory(
            str(PIPELINE), [str(PIPELINE / "config.yaml")]
        )
        second = build_source_inventory(
            str(PIPELINE), [str(PIPELINE / "config.yaml")]
        )
        self.assertEqual(first, second)

        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "source_files.sha256.tsv"
            self.assertTrue(
                write_source_hash_inventory(first, str(PIPELINE), str(output))
            )
            first_mtime = output.stat().st_mtime_ns
            time.sleep(0.02)
            self.assertFalse(
                write_source_hash_inventory(second, str(PIPELINE), str(output))
            )
            self.assertEqual(output.stat().st_mtime_ns, first_mtime)

            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(set(rows[0]), {"category", "sha256", "path"})
            paths = [row["path"] for row in rows]
            self.assertEqual(paths, sorted(paths))

    def test_missing_expected_source_fails_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            missing = root / "missing.py"
            with self.assertRaisesRegex(
                SourceInventoryError, "expected pipeline source is missing"
            ):
                write_source_hash_inventory(
                    [("helpers", str(missing))],
                    str(PIPELINE),
                    str(root / "hashes.tsv"),
                )


if __name__ == "__main__":
    unittest.main()

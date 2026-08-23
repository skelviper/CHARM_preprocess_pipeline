#!/usr/bin/env python3

import subprocess
from pathlib import Path
import unittest


PIPELINE = Path(__file__).resolve().parents[1]


class ReleaseDocumentationTests(unittest.TestCase):
    def test_pipeline_tree_excludes_external_runtime_and_run_outputs(self):
        self.assertFalse((PIPELINE / "vendor").exists())
        self.assertFalse((PIPELINE / "docs").exists())
        self.assertFalse((PIPELINE / "profiles").exists())
        self.assertFalse((PIPELINE / "slurm_log").exists())
        self.assertFalse((PIPELINE / "BASELINE_PROVENANCE.md").exists())
        self.assertFalse((PIPELINE / "MIGRATION_PROVENANCE.md").exists())
        self.assertFalse((PIPELINE / "SOURCE_MANIFEST.sha256").exists())
        self.assertTrue((PIPELINE.parent / "CHARMtools" / "UPSTREAM_MANIFEST.sha256").is_file())
        config = (PIPELINE / "config.yaml").read_text()
        self.assertIn("CHARMtools: ../CHARMtools", config)

    def test_launcher_help_exposes_only_per_cell_mapping(self):
        result = subprocess.run(
            ["bash", str(PIPELINE / "runCHARM.sh"), "--help"],
            cwd=str(PIPELINE),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("one BWA job per cell", result.stdout)
        self.assertIn("Snakemake's --cluster interface", result.stdout)
        self.assertIn("remain NOT_RUN", result.stdout)
        self.assertNotIn("mapping.mode", result.stdout)
        self.assertNotIn("shard", result.stdout.lower())

    def test_runnable_source_has_no_sharding_option(self):
        runtime_paths = [
            PIPELINE / "config.yaml",
            PIPELINE / "CHARM.smk",
            PIPELINE / "runCHARM.sh",
        ]
        runtime_paths.extend(sorted((PIPELINE / "rules").glob("*.rules")))
        runtime_paths.extend(sorted((PIPELINE / "CHARM_scripts").glob("*.py")))
        forbidden = (
            "mapping.mode",
            "mode: sharded",
            "--shard",
            "sharded_mapping",
            "shard_count",
        )
        for path in runtime_paths:
            text = path.read_text()
            for token in forbidden:
                self.assertNotIn(token, text, "{} references {}".format(path, token))

    def test_pipeline_uses_one_existing_conda_environment(self):
        runtime_paths = [PIPELINE / "CHARM.smk"]
        runtime_paths.extend(sorted((PIPELINE / "rules").glob("*.rules")))
        for path in runtime_paths:
            self.assertNotIn("conda:", path.read_text(), str(path))

        self.assertNotIn("--use-conda", (PIPELINE / "runCHARM.sh").read_text())
        self.assertEqual(
            [path.name for path in (PIPELINE / "envs").glob("*.yml")],
            ["charm.yml"],
        )

    def test_readme_describes_only_repository_purpose_and_io_layout(self):
        readme = (PIPELINE / "README.md").read_text()

        self.assertIn("## Purpose", readme)
        self.assertIn("## Input Layout", readme)
        self.assertIn("## Output Layout", readme)
        self.assertIn("<sample>_1.fq.gz", readme)
        self.assertIn("<sample>_R1.fq.gz", readme)
        self.assertIn("processed/", readme)
        self.assertIn("result/", readme)
        self.assertIn("qc/", readme)
        self.assertIn("tmp/", readme)
        self.assertNotIn("r1_all", readme)
        self.assertNotIn("r1_compatible", readme)
        self.assertNotIn("r1r2_concordant", readme)
        self.assertNotIn("Release candidate status", readme)
        self.assertNotIn("NOT_RUN", readme)

if __name__ == "__main__":
    unittest.main()

import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PIPELINE_DIR = Path(__file__).resolve().parents[1]
SCRIPT_DIR = PIPELINE_DIR / "CHARM_scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import input_contract


SNAKEMAKE = shutil.which("snakemake")


def write_fastq(path, read_id):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(str(path), "wt") as handle:
        handle.write("@{}\nACGT\n+\nIIII\n".format(read_id))


@unittest.skipUnless(SNAKEMAKE, "snakemake is not available in this environment")
class InputContractDagTest(unittest.TestCase):
    def test_dag_requires_a_current_launcher_generated_contract(self):
        with tempfile.TemporaryDirectory(prefix="charm-input-contract-dag-") as temporary:
            work_dir = Path(temporary)
            sample_dir = work_dir / "Rawdata" / "Alpha"
            write_fastq(sample_dir / "Alpha_1.fq.gz", "read_r1")
            write_fastq(sample_dir / "Alpha_2.fq.gz", "read_r2")
            target = "processed/Alpha/DNA/Alpha.dna.R1.fq.gz"
            command = [
                SNAKEMAKE,
                "--dry-run",
                "--cores",
                "1",
                "--snakefile",
                str(PIPELINE_DIR / "CHARM.smk"),
                target,
                "--config",
                "work_dir={}".format(work_dir),
            ]
            environment = os.environ.copy()
            environment["XDG_CACHE_HOME"] = str(work_dir / ".cache")

            missing = subprocess.run(
                command,
                cwd=str(work_dir),
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
            )
            self.assertNotEqual(missing.returncode, 0, missing.stdout)
            self.assertIn("invalid or missing Rawdata input contract", missing.stdout)

            input_contract.create_contract(str(work_dir))
            valid = subprocess.run(
                command,
                cwd=str(work_dir),
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
            )
            self.assertEqual(valid.returncode, 0, valid.stdout)
            self.assertIn("rule split", valid.stdout)

            write_fastq(sample_dir / "Alpha_1.fq.gz", "changed_r1")
            drifted = subprocess.run(
                command,
                cwd=str(work_dir),
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
            )
            self.assertNotEqual(drifted.returncode, 0, drifted.stdout)
            self.assertIn("Rawdata differs from the frozen input contract", drifted.stdout)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3

import copy
import gzip
import io
import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "CHARM_scripts"
import sys

sys.path.insert(0, str(SCRIPT_DIR))

import generate_stat_contract as STATS
import input_contract


def write_fastq(path, read_id, records=1):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(str(path), "wt") as handle:
        for index in range(records):
            handle.write("@{}_{}\nACGT\n+\nIIII\n".format(read_id, index))


def write_pairs(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(str(path), "wt") as handle:
        handle.write("#pairs format v1.0\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")


class GenerateStatContractTest(unittest.TestCase):
    def make_contract_workspace(self, root):
        backing = root / "backing"
        numeric = root / "Rawdata" / "Alpha_one"
        illumina = root / "Rawdata" / "Beta.two"
        write_fastq(backing / "alpha_r1.fq.gz", "alpha_r1", records=2)
        write_fastq(backing / "alpha_r2.fq.gz", "alpha_r2")
        numeric.mkdir(parents=True)
        os.symlink(
            os.path.relpath(backing / "alpha_r1.fq.gz", numeric),
            numeric / "Alpha_one_1.fq.gz",
        )
        os.symlink(
            os.path.relpath(backing / "alpha_r2.fq.gz", numeric),
            numeric / "Alpha_one_2.fq.gz",
        )
        write_fastq(illumina / "Beta.two_R1.fq.gz", "beta_r1", records=3)
        write_fastq(illumina / "Beta.two_R2.fq.gz", "beta_r2")
        contract = input_contract.create_contract(str(root))
        return contract, root / "qc" / "input_contract" / "current.json"

    def test_manifest_and_raw_stats_use_both_frozen_naming_schemes(self):
        with tempfile.TemporaryDirectory() as temporary_name:
            root = Path(temporary_name)
            contract, contract_path = self.make_contract_workspace(root)
            loaded = STATS.load_frozen_contract(contract_path, root)
            manifest = root / "manifest.tsv"
            raw_stats = root / "raw.fq.stat"
            STATS.write_manifest(loaded, manifest)
            STATS.write_raw_fastq_stats(loaded, raw_stats, workers=2)

            text = manifest.read_text()
            self.assertIn("Alpha_one\t", text)
            self.assertEqual("sample\tsafe_code", text.splitlines()[0])
            self.assertEqual(
                [
                    "./Rawdata/Alpha_one/Alpha_one_1.fq.gz\t8",
                    "./Rawdata/Beta.two/Beta.two_R1.fq.gz\t12",
                ],
                raw_stats.read_text().splitlines(),
            )
            self.assertEqual(2, contract["sample_count"])

    def test_worker_keeps_frozen_cohort_and_resolved_symlink_target(self):
        with tempfile.TemporaryDirectory() as temporary_name:
            root = Path(temporary_name)
            _, contract_path = self.make_contract_workspace(root)
            frozen = STATS.load_frozen_contract(contract_path, root)

            drift = root / "Rawdata" / "Gamma"
            write_fastq(drift / "Gamma_1.fq.gz", "gamma_r1", records=7)
            write_fastq(drift / "Gamma_2.fq.gz", "gamma_r2")
            replacement = root / "backing" / "replacement.fq.gz"
            write_fastq(replacement, "replacement", records=9)
            alpha_link = root / "Rawdata" / "Alpha_one" / "Alpha_one_1.fq.gz"
            alpha_link.unlink()
            os.symlink(os.path.relpath(replacement, alpha_link.parent), alpha_link)

            output = root / "raw.fq.stat"
            STATS.write_raw_fastq_stats(frozen, output, workers=2)
            self.assertEqual(2, len(output.read_text().splitlines()))
            self.assertIn("Alpha_one_1.fq.gz\t8", output.read_text())
            self.assertNotIn("Gamma", output.read_text())
            with self.assertRaises(input_contract.InputContractError):
                input_contract.load_and_validate_contract(str(root))

    def test_rna_qname_uses_penultimate_safe_code_and_authoritative_name(self):
        with tempfile.TemporaryDirectory() as temporary_name:
            root = Path(temporary_name)
            contract, _ = self.make_contract_workspace(root)
            mapping = STATS.code_to_sample(contract)
            code = next(
                sample["safe_code"]
                for sample in contract["samples"]
                if sample["sample_name"] == "Alpha_one"
            )
            qname = "instrument_lane_read_with_underscores_{}_ACGTACGT".format(code)
            sam = (
                "{}\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tXS:Z:Assigned\n"
                "{}\t256\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
            ).format(qname, qname)
            output = root / "rna.stat"
            STATS.write_rna_alignment_stats(contract, io.StringIO(sam), output)
            self.assertEqual("Alpha_one\t2\t1\n", output.read_text())
            self.assertEqual("Alpha_one", STATS.sample_from_qname(qname, mapping))

            for invalid in (
                "too_short",
                "prefix_C0000000000000000_BADUMI",
                "prefix_Cffffffffffffffff_ACGTACGT",
            ):
                with self.assertRaises(STATS.StatisticsContractError):
                    STATS.sample_from_qname(invalid, mapping)

    def test_duplicate_contract_codes_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary_name:
            root = Path(temporary_name)
            contract, contract_path = self.make_contract_workspace(root)
            tampered = copy.deepcopy(contract)
            tampered["samples"][1]["safe_code"] = tampered["samples"][0]["safe_code"]
            tampered["contract_sha256"] = input_contract._payload_sha256(
                input_contract._contract_core(tampered)
            )
            contract_path.write_text(json.dumps(tampered))
            with self.assertRaises(input_contract.InputContractError):
                STATS.load_frozen_contract(contract_path, root)

    def test_c123_and_rmsd_stats_have_separate_complete_contracts(self):
        with tempfile.TemporaryDirectory() as temporary_name:
            root = Path(temporary_name)
            contract, _ = self.make_contract_workspace(root)
            for sample in ("Alpha_one", "Beta.two"):
                write_pairs(
                    root / "result" / "cleaned_pairs" / "c123" / "{}.pairs.gz".format(sample),
                    [
                        ["r1", "chr1", "10", "chr1", "20"],
                        ["r2", "chr1", "30", "chr2", "40"],
                    ],
                )
                rmsd = (
                    root
                    / "result"
                    / "3d_info"
                    / sample
                    / "{}.20k.align.rms.info".format(sample)
                )
                rmsd.parent.mkdir(parents=True, exist_ok=True)
                rmsd.write_text("[M::__main__] top3 RMS RMSD: 1.25\n")

            pairs_output = root / "pairs.stat"
            inter_output = root / "inter.stat"
            rmsd_output = root / "rmsd.stat"
            STATS.write_cleaned_pair_stats(
                contract,
                root,
                pairs_output,
                inter_output,
                workers=2,
            )
            self.assertTrue(all(line.endswith("\t2") for line in pairs_output.read_text().splitlines()))
            self.assertTrue(all(line.endswith("\t1") for line in inter_output.read_text().splitlines()))
            STATS.write_rmsd_stats(contract, root, ["20k"], rmsd_output)
            self.assertEqual(2, len(rmsd_output.read_text().splitlines()))

            (root / "result" / "3d_info" / "Beta.two" / "Beta.two.20k.align.rms.info").unlink()
            with self.assertRaisesRegex(
                STATS.StatisticsContractError, "missing required output"
            ):
                STATS.write_rmsd_stats(contract, root, ["20k"], rmsd_output)

    def test_shell_has_strict_contract_bound_discovery(self):
        text = (SCRIPT_DIR / "generateStat.sh").read_text()
        self.assertIn("set -euo pipefail", text)
        self.assertIn("generate_stat_contract.py", text)
        self.assertIn('"qc/stat/$output_name"', text)
        self.assertIn("charm|hires", text)
        for forbidden in ("ls Rawdata", "find -L ./Rawdata", "*1.fq.gz"):
            self.assertNotIn(forbidden, text)

    def test_notebook_reads_qc_local_stat_and_writes_metadata(self):
        text = (SCRIPT_DIR / "stat.ipynb").read_text()
        self.assertNotIn("../stat/", text)
        self.assertNotIn("if_charm", text)
        self.assertNotIn("/mnt/ssd/zliu/run_charm", text)
        self.assertIn("provenance/effective_config.json", text)
        self.assertIn('experiment_type == \\\"charm\\\"', text)
        self.assertIn("if (config$if_structure)", text)
        self.assertIn("stat/raw.fq.stat", text)
        self.assertNotIn('strip_after = \\\"_\\\"', text)
        self.assertIn("input_contract/discovered_cells.tsv", text)
        self.assertIn("metadata cell set does not match the frozen input contract", text)
        self.assertIn("experiment_type = experiment_type", text)
        self.assertIn("RNA_primary_output_type = primary_rna_output_type", text)
        self.assertIn("primary_rna_output_type <- config$rna_primary_output_type", text)
        self.assertIn("UMIs_gene = colSums(rna_gene_matrix)", text)
        self.assertIn("UMIs_exon = colSums(rna_exon_matrix)", text)
        self.assertIn("stat/pairs.c123.stat", text)
        self.assertIn("pairsValidRatio = ifelse(raw_pairs > 0", text)
        self.assertIn("interPairsRatio = ifelse(pairs_clean3 > 0", text)
        self.assertNotIn("feature_stat_modes", text)
        self.assertNotIn("rna_compatibility_path", text)
        self.assertNotIn("UMIs_gene_r1_all", text)
        self.assertIn("canonical_ensembl_levels", text)
        self.assertNotIn("seqlevelsStyle(annotations)", text)
        self.assertNotIn("%>% slice(", text)
        self.assertIn("%>% dplyr::slice(", text)
        self.assertIn('write_tsv(\\\"metadata_raw.tsv\\\")', text)

    def test_final_audit_keeps_mode_metrics_outside_primary_metadata(self):
        text = (SCRIPT_DIR / "audit_complete_run.py").read_text()
        self.assertIn('"rna_mode_tables"', text)
        self.assertIn("rna.output_modes.per_cell.tsv", text)
        self.assertIn("rna.output_modes.summary.tsv", text)
        self.assertIn("counts.{}.total.format.tsv", text)
        self.assertIn("non-primary RNA columns leaked into metadata", text)
        self.assertIn('umi_column = "UMIs_{}".format(feature)', text)
        self.assertIn('"pairs_clean3"', text)
        self.assertIn("metadata/c123 count mismatch", text)
        self.assertIn("metadata/c123 inter-ratio mismatch", text)

    def test_statistics_shell_publishes_no_disabled_structure_placeholders(self):
        with tempfile.TemporaryDirectory() as temporary_name:
            root = Path(temporary_name)
            contract, contract_path = self.make_contract_workspace(root)
            for sample_entry in contract["samples"]:
                sample = sample_entry["sample_name"]
                write_fastq(
                    root / "processed" / sample / "RNA" / "{}.rna.clean.R1.fq.gz".format(sample),
                    "rna_r1",
                )
                write_fastq(
                    root / "processed" / sample / "RNA" / "{}.rna.clean.R2.fq.gz".format(sample),
                    "rna_r2",
                )
                write_fastq(
                    root / "processed" / sample / "DNA" / "{}.dna.clean.R1.fq.gz".format(sample),
                    "dna_r1",
                )
                two_d = root / "processed" / sample / "2d_info"
                write_pairs(
                    two_d / "raw.pairs.gz",
                    [["r1", "chr1", "1", "chr1", "2"]],
                )
                write_pairs(
                    two_d / "contacts.pairs.gz",
                    [["r1", "chr1", "1", "chr1", "2"]],
                )
                (two_d / "{}.yperx.txt".format(sample)).write_text("0.01\n")
                for clean in ("c1", "c12", "c123"):
                    write_pairs(
                        root
                        / "result"
                        / "cleaned_pairs"
                        / clean
                        / "{}.pairs.gz".format(sample),
                        [["r1", "chr1", "1", "chr1", "2"]],
                    )
                for split in ("atac", "ct"):
                    bam = root / "processed" / sample / split / "{}.sort.bam".format(sample)
                    bam.parent.mkdir(parents=True, exist_ok=True)
                    bam.write_bytes(b"fixture")
                    fragment_log = (
                        root
                        / "processed"
                        / "{}_all".format(split)
                        / "{}.{}.frag.log".format(sample, split)
                    )
                    fragment_log.parent.mkdir(parents=True, exist_ok=True)
                    fragment_log.write_text("Combined duplication rate is 12.50%\n")

            fake_bin = root / "fake-bin"
            fake_bin.mkdir()
            rna_count_bam = (
                root
                / "processed"
                / "RNA_all"
                / "feature_r1_all_gene_total"
                / "samsort.bam"
            )
            rna_count_bam.parent.mkdir(parents=True)
            rna_count_bam.write_bytes(b"fixture")
            safe_codes = [sample["safe_code"] for sample in contract["samples"]]
            fake_samtools = fake_bin / "samtools"
            fake_samtools.write_text(
                "#!/usr/bin/env bash\n"
                "set -euo pipefail\n"
                "case $1 in\n"
                "  flagstat) printf '10 + 0 in total (QC-passed reads + QC-failed reads)\\n' ;;\n"
                "  view)\n"
                "    printf 'read_{}_ACGTACGT\\t0\\tchr1\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\tIIII\\tXS:Z:Assigned\\n'\n"
                "    printf 'read_{}_ACGTACGT\\t0\\tchr1\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\tIIII\\tXS:Z:Assigned\\n'\n"
                "    ;;\n"
                "  *) exit 2 ;;\n"
                "esac\n".format(safe_codes[0], safe_codes[1])
            )
            fake_samtools.chmod(0o755)
            environment = os.environ.copy()
            environment["PATH"] = str(fake_bin) + os.pathsep + environment["PATH"]
            environment["CHARM_STAT_THREADS"] = "2"
            command = [
                "bash",
                str(SCRIPT_DIR / "generateStat.sh"),
                str(root),
                str(contract_path),
                "False",
                "charm",
                "20k,50k,200k,1m",
                str(rna_count_bam),
            ]
            disabled = subprocess.run(
                command,
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            self.assertEqual(0, disabled.returncode, disabled.stdout)
            self.assertTrue((root / "qc" / "stat" / "pairs.c123.stat").is_file())
            self.assertTrue((root / "qc" / "stat" / "inter.pairs.c123.stat").is_file())
            self.assertFalse((root / "qc" / "stat" / "rmsd.info").exists())
            self.assertEqual(
                {"Alpha_one,10", "Beta.two,10"},
                set((root / "qc" / "stat" / "atac.read.stat").read_text().splitlines()),
            )
            self.assertEqual(
                {"Alpha_one\t1\t1", "Beta.two\t1\t1"},
                set(
                    (root / "qc" / "stat" / "rna.reads_per_cell.stat")
                    .read_text()
                    .splitlines()
                ),
            )

            shutil.rmtree(root / "qc" / "stat")
            command[5] = "hires"
            hires = subprocess.run(
                command,
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            self.assertEqual(0, hires.returncode, hires.stdout)
            self.assertTrue((root / "qc" / "stat" / "raw.fq.stat").is_file())
            self.assertFalse((root / "qc" / "stat" / "atac.read.stat").exists())
            self.assertFalse((root / "qc" / "stat" / "ct.read.stat").exists())

            shutil.rmtree(root / "qc" / "stat")
            command[4] = "True"
            enabled_missing = subprocess.run(
                command,
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            self.assertNotEqual(0, enabled_missing.returncode, enabled_missing.stdout)
            self.assertIn("missing required output", enabled_missing.stdout)
            self.assertFalse((root / "qc" / "stat").exists())


if __name__ == "__main__":
    unittest.main()

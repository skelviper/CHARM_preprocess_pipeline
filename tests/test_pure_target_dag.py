#!/usr/bin/env python3

import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

import yaml


PIPELINE = Path(__file__).resolve().parents[1]
SNAKEFILE = PIPELINE / "CHARM.smk"
INPUT_CONTRACT = PIPELINE / "CHARM_scripts" / "input_contract.py"
BWA_INDEX_SUFFIXES = (
    "",
    ".0123",
    ".amb",
    ".ann",
    ".bwt.2bit.64",
    ".pac",
    ".fai",
)
STAR_INDEX_COMPONENTS = (
    "Genome",
    "Log.out",
    "SA",
    "SAindex",
    "chrLength.txt",
    "chrName.txt",
    "chrNameLength.txt",
    "chrStart.txt",
    "exonGeTrInfo.tab",
    "exonInfo.tab",
    "geneInfo.tab",
    "genomeParameters.txt",
    "sjdbInfo.txt",
    "sjdbList.fromGTF.out.tab",
    "sjdbList.out.tab",
    "transcriptInfo.tab",
)


class PureTargetDagTests(unittest.TestCase):
    def setUp(self):
        if shutil.which("snakemake") is None:
            self.skipTest("snakemake is not available in PATH")
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.workspace = self.root / "workspace"
        self.workspace.mkdir()
        self._write_raw_fastqs("S1")
        self.config_path = self.root / "config.yaml"
        self._write_config()
        subprocess.run(
            [
                sys.executable,
                str(INPUT_CONTRACT),
                "create",
                "--work-dir",
                str(self.workspace),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    def tearDown(self):
        self.temporary.cleanup()

    def _write_raw_fastqs(self, sample):
        sample_dir = self.workspace / "Rawdata" / sample
        sample_dir.mkdir(parents=True)
        records = {
            1: "@read1/1\nACGTACGT\n+\nIIIIIIII\n",
            2: "@read1/2\nTGCATGCA\n+\nIIIIIIII\n",
        }
        for read, record in records.items():
            with gzip.open(str(sample_dir / "{}_{}.fq.gz".format(sample, read)), "wt") as handle:
                handle.write(record)

    def _write_config(self):
        with open(str(PIPELINE / "config.yaml")) as handle:
            config = yaml.safe_load(handle)

        reference_root = self.root / "reference"
        star_index = reference_root / "star"
        star_index.mkdir(parents=True)
        bwa_reference = reference_root / "genome.fa"
        annotation = reference_root / "annotation.gtf"
        snp = reference_root / "snp.tsv"
        par = reference_root / "par.bed"
        rna_snp = reference_root / "rna_snp.tsv.gz"
        for suffix in BWA_INDEX_SUFFIXES:
            Path(str(bwa_reference) + suffix).write_text(
                ">chr1\nACGT\n" if not suffix else "fixture\n"
            )
        for component in STAR_INDEX_COMPONENTS:
            (star_index / component).write_text("fixture\n")
        annotation.write_text(
            'chr1\ttest\tgene\t1\t4\t.\t+\t.\tgene_id "g1"; gene_name "g1";\n'
        )
        snp.write_text("chr1\t1\tA\tG\n")
        par.write_text("chr1\t0\t1\n")
        with gzip.open(str(rna_snp), "wt") as handle:
            handle.write("chr1\t1\tA\tG\n")

        config["work_dir"] = str(self.workspace)
        config["if_structure"] = False
        config["experiment_type"] = "charm"
        config["if_RNA_snp_split"] = False
        # A stale external setting must not restore the former no-op behavior.
        config["relaxed_rule_all"] = True
        config["refs"][config["ref_genome"]] = {
            "bwa_mem2_index": str(bwa_reference),
            "star_index": str(star_index),
            "annotations": str(annotation),
            "snp": str(snp),
            "par": str(par),
            "RNAsnp": str(rna_snp),
        }
        config["cuttag_fragments"]["temp_root"] = str(self.workspace / "tmp")
        with open(str(self.config_path), "w") as handle:
            yaml.safe_dump(config, handle, sort_keys=False)

    def _update_config(self, **updates):
        with open(str(self.config_path)) as handle:
            config = yaml.safe_load(handle)
        config.update(updates)
        with open(str(self.config_path), "w") as handle:
            yaml.safe_dump(config, handle, sort_keys=False)

    def _snakemake(self, *targets, dry_run=False, extra=None):
        command = [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "--configfile",
            str(self.config_path),
            "--cores",
            "1",
            "--nolock",
        ]
        if dry_run:
            command.append("--dry-run")
        if extra:
            command.extend(extra)
        command.extend(targets)
        return subprocess.run(
            command,
            cwd=str(PIPELINE),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

    def test_module_targets_dry_run_independently(self):
        dry_runs = {}
        for target in ("rna", "contacts2d", "fragments", "provenance", "qc", "all"):
            result = self._snakemake(target, dry_run=True)
            self.assertEqual(result.returncode, 0, "{}\n{}".format(target, result.stdout))
            self.assertNotIn("Nothing to be done", result.stdout, target)
            dry_runs[target] = result.stdout

        fragments = dry_runs["fragments"]
        self.assertEqual(fragments.count("rule multiSplit:"), 1)
        self.assertIn("processed/S1/ct/S1.ct.R1.fq.gz", fragments)
        self.assertIn("processed/S1/atac/S1.atac.R1.fq.gz", fragments)

        contacts2d = dry_runs["contacts2d"]
        self.assertIn("rule clean_leg:", contacts2d)
        self.assertIn("rule clean_isolated:", contacts2d)
        self.assertIn("rule clean_splicing:", contacts2d)

        qc = dry_runs["qc"]
        self.assertIn("rule generate_statistics:", qc)
        self.assertIn("rule stage_qc_notebook:", qc)
        self.assertIn("rule generate_metadata:", qc)
        self.assertIn("rule audit_complete_run:", qc)
        self.assertIn("qc/stat/raw.fq.stat", qc)
        self.assertIn("qc/stat.ipynb", qc)
        self.assertIn("qc/stat.executed.ipynb", qc)
        self.assertIn("qc/metadata_raw.tsv", qc)
        self.assertIn("qc/COMPLETE_RUN_AUDIT.tsv", qc)

        overall = dry_runs["all"]
        self.assertIn("rule generate_statistics:", overall)
        self.assertIn("rule generate_metadata:", overall)
        self.assertIn("rule audit_complete_run:", overall)
        self.assertIn("qc/target_outputs.tsv", overall)
        self.assertNotIn("aggregate_qc_receipts", overall)

    def test_missing_outputs_cannot_use_relaxed_noop(self):
        result = self._snakemake("all", dry_run=True)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("rule all", result.stdout)
        self.assertNotIn("Nothing to be done", result.stdout)
        self.assertFalse((self.workspace / "qc" / "stat").exists())
        self.assertFalse((self.workspace / "qc" / "stat.ipynb").exists())
        self.assertTrue(
            (self.workspace / "qc" / "input_contract" / "current.json").exists()
        )
        self.assertFalse((self.workspace / "qc" / "target_outputs.tsv").exists())

    def test_structure_rules_and_target_are_absent_when_disabled(self):
        listed = self._snakemake(extra=["--list-target-rules"])
        self.assertEqual(listed.returncode, 0, listed.stdout)
        target_rules = set(listed.stdout.splitlines())
        self.assertNotIn("structure3d", target_rules)

        all_result = self._snakemake("all", dry_run=True)
        self.assertEqual(all_result.returncode, 0, all_result.stdout)
        for rule_name in (
            "hickit_3d",
            "compress_3dg",
            "hickit_clean3D",
            "compress_clean3d",
            "hickit_align3D",
        ):
            self.assertNotIn(rule_name, all_result.stdout)

        disabled = self._snakemake("structure3d", dry_run=True)
        self.assertNotEqual(disabled.returncode, 0, disabled.stdout)

        qc = self._snakemake("qc", dry_run=True)
        self.assertEqual(qc.returncode, 0, qc.stdout)
        for enabled_2d_stat in (
            "qc/stat/pairs.c123.stat",
            "qc/stat/inter.pairs.c123.stat",
        ):
            self.assertIn(enabled_2d_stat, qc.stdout)
        self.assertIn("rule clean_splicing:", qc.stdout)
        self.assertNotIn("qc/stat/rmsd.info", qc.stdout)

    def test_hires_omits_charm_fragments_and_charm_statistics(self):
        self._update_config(experiment_type="hires")

        qc = self._snakemake("qc", dry_run=True)
        self.assertEqual(qc.returncode, 0, qc.stdout)
        self.assertIn("qc/stat/raw.fq.stat", qc.stdout)
        for charm_only in (
            "qc/stat/atac.read.stat",
            "qc/stat/ct.read.stat",
            "qc/stat/atac.dedup_rate.stat",
            "qc/stat/ct.dedup_rate.stat",
            "rule build_cuttag_phased_snp_index:",
            "rule bam2bed:",
        ):
            self.assertNotIn(charm_only, qc.stdout)

        contract_result = self._snakemake("target_contract")
        self.assertEqual(contract_result.returncode, 0, contract_result.stdout)
        contract = (self.workspace / "qc" / "target_outputs.tsv").read_text()
        self.assertIn("fragments\t0\tNA", contract)
        self.assertNotIn("qc/stat/atac.read.stat", contract)
        self.assertNotIn("qc/stat/ct.read.stat", contract)

    def test_invalid_experiment_type_fails_before_dag_construction(self):
        self._update_config(experiment_type="other")
        result = self._snakemake("all", dry_run=True)
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("experiment_type must be exactly 'charm' or 'hires'", result.stdout)

    def test_rna_target_expands_all_configured_output_types(self):
        result = self._snakemake("rna", dry_run=True)
        self.assertEqual(result.returncode, 0, result.stdout)
        for mode in ("r1_all", "r1_compatible", "r1r2_concordant"):
            self.assertIn(
                "result/RNA_Res/{}/counts.gene.total.tsv".format(mode),
                result.stdout,
            )
        self.assertIn("rule star_mapping_r2_validator:", result.stdout)
        rna_rules = (PIPELINE / "rules" / "CHARM_RNA.rules").read_text()
        self.assertIn("--poly-a --minimum-length 0:1", rna_rules)
        self.assertIn("contract_file=$(realpath {input.input_contract})", rna_rules)
        self.assertNotIn('{WORK_DIR}/{input.input_contract}', rna_rules)
        self.assertIn("rule filter_rna_gene_compatibility:", result.stdout)
        self.assertIn("rule summarize_rna_output_modes:", result.stdout)

    def test_rna_heavy_jobs_expose_serial_resource_slots(self):
        rna_rules = (PIPELINE / "rules" / "CHARM_RNA.rules").read_text()
        self.assertEqual(rna_rules.count("star_slots=1"), 2)
        self.assertEqual(rna_rules.count("count_slots=1"), 2)

    def test_r1_all_only_skips_r2_validator_branch(self):
        self._update_config(rna_output_types=["r1_all"])
        result = self._snakemake("rna", dry_run=True)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("result/RNA_Res/r1_all/counts.gene.total.tsv", result.stdout)
        self.assertNotIn("r1_compatible/counts", result.stdout)
        self.assertNotIn("r1r2_concordant/counts", result.stdout)
        self.assertNotIn("rule star_mapping_r2_validator:", result.stdout)
        self.assertNotIn("rule filter_rna_gene_compatibility:", result.stdout)

    def test_invalid_rna_output_type_contract_fails_dag_construction(self):
        for values, message in (
            ([], "must be a non-empty list"),
            (["r1_all", "r1_all"], "must not contain duplicate"),
            (["r1_all", "unknown"], "invalid rna_output_types"),
        ):
            self._update_config(rna_output_types=values)
            result = self._snakemake("rna", dry_run=True)
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn(message, result.stdout)

    def test_primary_rna_output_type_must_be_selected(self):
        self._update_config(
            rna_output_types=["r1_all"],
            rna_primary_output_type="r1_compatible",
        )
        result = self._snakemake("all", dry_run=True)
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn(
            "rna_primary_output_type must be one of the selected rna_output_types",
            result.stdout,
        )

    def test_structure_target_closes_when_enabled(self):
        self._update_config(if_structure=True)

        listed = self._snakemake(extra=["--list-target-rules"])
        self.assertEqual(listed.returncode, 0, listed.stdout)
        self.assertIn("structure3d", set(listed.stdout.splitlines()))

        structure = self._snakemake("structure3d", dry_run=True)
        self.assertEqual(structure.returncode, 0, structure.stdout)
        for rule_name in (
            "clean_splicing",
            "hickit_3d",
            "hickit_clean3D",
            "compress_clean3d",
            "hickit_align3D",
        ):
            self.assertIn(rule_name, structure.stdout)

        qc = self._snakemake("qc", dry_run=True)
        self.assertEqual(qc.returncode, 0, qc.stdout)
        self.assertIn("rule generate_statistics:", qc.stdout)
        self.assertIn("rule clean_splicing:", qc.stdout)
        for enabled_stat in (
            "qc/stat/pairs.c123.stat",
            "qc/stat/inter.pairs.c123.stat",
            "qc/stat/rmsd.info",
        ):
            self.assertIn(enabled_stat, qc.stdout)

    def test_rule_all_is_declarative(self):
        text = SNAKEFILE.read_text()
        all_start = text.index("rule all:")
        all_end = text.index("\n\nrule rna:", all_start)
        all_block = text[all_start:all_end]
        for forbidden in ("shell:", "run:", "mkdir", "cp ", "generateStat"):
            self.assertNotIn(forbidden, all_block)
        self.assertNotIn("RELAXED_RULE_ALL", text)
        self.assertNotIn("_rule_all_inputs", text)
        self.assertNotIn("relaxed_rule_all:", (PIPELINE / "config.yaml").read_text())
        self.assertNotIn("stat.ipynb", (PIPELINE / "runCHARM.sh").read_text())

    def test_target_contract_is_machine_readable(self):
        result = self._snakemake("target_contract")
        self.assertEqual(result.returncode, 0, result.stdout)
        contract = self.workspace / "qc" / "target_outputs.tsv"
        with open(str(contract), newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(set(rows[0]), {"target", "enabled", "output"})
        by_target = {}
        for row in rows:
            by_target.setdefault(row["target"], []).append(row)
        self.assertEqual(
            set(by_target),
            {"rna", "contacts2d", "fragments", "structure3d", "provenance", "qc", "all"},
        )
        self.assertEqual(
            by_target["structure3d"],
            [{"target": "structure3d", "enabled": "0", "output": "NA"}],
        )
        qc_outputs = {row["output"] for row in by_target["qc"]}
        self.assertIn("qc/stat/pairs.c123.stat", qc_outputs)
        self.assertIn("qc/stat/inter.pairs.c123.stat", qc_outputs)
        self.assertNotIn("qc/stat/rmsd.info", qc_outputs)
        contacts_outputs = {row["output"] for row in by_target["contacts2d"]}
        self.assertIn("result/cleaned_pairs/c123/S1.pairs.gz", contacts_outputs)
        self.assertIn(
            "result/fragments/atac.fragments.bgz.tbi",
            {row["output"] for row in by_target["fragments"]},
        )
        self.assertIn(
            "qc/provenance/effective_config.json",
            {row["output"] for row in by_target["all"]},
        )
        self.assertIn("qc/stat/raw.fq.stat", qc_outputs)
        self.assertIn("qc/stat.executed.ipynb", qc_outputs)
        self.assertIn("qc/metadata_raw.tsv", qc_outputs)
        self.assertIn("qc/COMPLETE_RUN_AUDIT.tsv", qc_outputs)
        rna_outputs = {row["output"] for row in by_target["rna"]}
        self.assertIn(
            "result/RNA_Res/r1_all/counts.gene.total.tsv", rna_outputs
        )
        self.assertIn(
            "result/RNA_Res/r1_compatible/counts.gene.total.tsv", rna_outputs
        )
        self.assertIn(
            "result/RNA_Res/r1r2_concordant/counts.gene.total.tsv", rna_outputs
        )
        self.assertIn("qc/stat/rna.output_modes.summary.tsv", qc_outputs)

    def test_target_contract_tracks_config_changes(self):
        initial = self._snakemake("target_contract")
        self.assertEqual(initial.returncode, 0, initial.stdout)

        self._update_config(if_structure=True)
        changed = self._snakemake("target_contract", dry_run=True)
        self.assertEqual(changed.returncode, 0, changed.stdout)
        self.assertIn("rule target_contract", changed.stdout)

        updated = self._snakemake("target_contract")
        self.assertEqual(updated.returncode, 0, updated.stdout)
        contract = (self.workspace / "qc" / "target_outputs.tsv").read_text()
        self.assertIn("result/3d_info/S1/S1.20k.align.rms.info", contract)
        self.assertIn("qc\t1\tqc/stat/pairs.c123.stat", contract)
        self.assertIn("qc\t1\tqc/stat/inter.pairs.c123.stat", contract)
        self.assertIn("qc\t1\tqc/stat/rmsd.info", contract)

    def test_provenance_tracks_effective_config_and_source_hashes(self):
        initial = self._snakemake("provenance")
        self.assertEqual(initial.returncode, 0, initial.stdout)

        provenance_dir = self.workspace / "qc" / "provenance"
        effective_config = json.loads(
            (provenance_dir / "effective_config.json").read_text()
        )
        self.assertEqual(effective_config["work_dir"], str(self.workspace))

        with open(str(provenance_dir / "source_files.sha256.tsv"), newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        source_hashes = {row["path"]: row["sha256"] for row in rows}
        self.assertEqual(
            source_hashes["CHARM.smk"],
            hashlib.sha256(SNAKEFILE.read_bytes()).hexdigest(),
        )
        classifier = PIPELINE / "CHARM_scripts" / "classify_2d_bam.py"
        self.assertEqual(
            source_hashes["CHARM_scripts/classify_2d_bam.py"],
            hashlib.sha256(classifier.read_bytes()).hexdigest(),
        )

        self._update_config(yperx_threshold=0.123)
        changed = self._snakemake("provenance", dry_run=True)
        self.assertEqual(changed.returncode, 0, changed.stdout)
        self.assertIn("rule pipeline_provenance", changed.stdout)

        updated = self._snakemake("provenance")
        self.assertEqual(updated.returncode, 0, updated.stdout)
        effective_config = json.loads(
            (provenance_dir / "effective_config.json").read_text()
        )
        self.assertEqual(effective_config["yperx_threshold"], 0.123)

    def test_completed_all_target_dry_runs_to_zero_jobs(self):
        summary = self._snakemake("all", extra=["--summary"])
        self.assertEqual(summary.returncode, 0, summary.stdout)
        output_paths = []
        for line in summary.stdout.splitlines():
            if not line or line.startswith("output_file"):
                continue
            fields = line.split("\t")
            if len(fields) >= 2 and fields[0] != "-":
                output_paths.append(fields[0])
        self.assertGreater(len(output_paths), 20, summary.stdout)
        completed_ns = time.time_ns()
        for output_path in output_paths:
            path = Path(output_path)
            if not path.is_absolute():
                path = self.workspace / path
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
            # The summary is not guaranteed to be topologically ordered. Give
            # every simulated workflow output one completion timestamp so no
            # generated prerequisite appears newer than its consumer.
            os.utime(str(path), ns=(completed_ns, completed_ns))

        rerun = self._snakemake("all", dry_run=True)
        self.assertEqual(rerun.returncode, 0, rerun.stdout)
        self.assertIn("Nothing to be done", rerun.stdout)


if __name__ == "__main__":
    unittest.main()

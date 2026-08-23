#!/usr/bin/env python3

import gzip
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "2d"
AUDIT_PATH = ROOT / "tests" / "audit_2d_equivalence.py"
SPEC = importlib.util.spec_from_file_location("audit_2d_equivalence", str(AUDIT_PATH))
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)
INPUT_CONTRACT = ROOT / "CHARM_scripts" / "input_contract.py"
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


def run(command, **kwargs):
    return subprocess.run(command, check=True, text=True, **kwargs)


def build_dag_fixture(root):
    workspace = root / "workspace"
    sample_dir = workspace / "Rawdata" / "fixture"
    sample_dir.mkdir(parents=True)
    for mate in (1, 2):
        with gzip.open(str(sample_dir / "fixture_{}.fq.gz".format(mate)), "wt") as handle:
            handle.write("@read/{}\nACGT\n+\nIIII\n".format(mate))

    reference = root / "reference"
    star_index = reference / "star"
    star_index.mkdir(parents=True)
    bwa_reference = reference / "genome.fa"
    annotation = reference / "annotation.gtf"
    snp = reference / "snp.tsv"
    par = reference / "par.bed"
    rna_snp = reference / "rna_snp.tsv.gz"
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

    with open(str(ROOT / "config.yaml")) as handle:
        config = yaml.safe_load(handle)
    config["work_dir"] = str(workspace)
    config["if_structure"] = False
    config["refs"][config["ref_genome"]] = {
        "bwa_mem2_index": str(bwa_reference),
        "star_index": str(star_index),
        "annotations": str(annotation),
        "snp": str(snp),
        "par": str(par),
        "RNAsnp": str(rna_snp),
    }
    config_path = root / "config.yaml"
    with open(str(config_path), "w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)

    run(
        [
            sys.executable,
            str(INPUT_CONTRACT),
            "create",
            "--work-dir",
            str(workspace),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return workspace, config_path


class RuleStructureTests(unittest.TestCase):
    def test_direct_name_sort_and_declarative_clean_rules(self):
        rules = (ROOT / "rules" / "scHiC_2dprocess.rules").read_text()
        self.assertIn("classify_2d_bam.py", rules)
        self.assertNotIn("samtools view -q 30 {input.bam} chrX chrY", rules)
        self.assertNotIn("charm_bam_records", rules)
        self.assertEqual(rules.count("-q {params.contact_min_mapq}"), 4)
        self.assertIn('if [ "$sex_call" = "XY" ]', rules)
        self.assertIn("{log.sex_classification}", rules)
        self.assertNotIn("samtools view -h {input.bam}", rules)
        self.assertIn("samtools sort -n --no-PG -O SAM", rules)
        self.assertNotIn("samtools collate", rules)
        self.assertIn("rule clean_leg:", rules)
        self.assertIn("rule clean_isolated:", rules)
        self.assertIn("rule clean_splicing:", rules)

    @unittest.skipUnless(shutil.which("snakemake"), "snakemake is not on PATH")
    def test_c12_target_does_not_schedule_c123(self):
        with tempfile.TemporaryDirectory() as tmp:
            workspace, config_path = build_dag_fixture(Path(tmp))
            source = workspace / "processed" / "fixture" / "2d_info"
            source.mkdir(parents=True)
            with gzip.open(str(source / "contacts.pairs.gz"), "wt") as handle:
                handle.write("## pairs format v1.0\n")
                handle.write("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n")
                handle.write("q1\tchr1\t1000\tchr1\t50000\t+\t+\n")
            result = subprocess.run(
                [
                    "snakemake",
                    "--dry-run",
                    "--cores",
                    "1",
                    "--snakefile",
                    str(ROOT / "CHARM.smk"),
                    "--configfile",
                    str(config_path),
                    "--nolock",
                    "result/cleaned_pairs/c12/fixture.pairs.gz",
                ],
                cwd=str(ROOT),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertIn("clean_leg", result.stdout)
            self.assertIn("clean_isolated", result.stdout)
            self.assertNotIn("clean_splicing", result.stdout)
            self.assertNotIn("hickit_2d", result.stdout)


class ContactParityTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which("samtools"), "samtools is not on PATH")
    def test_baseline_and_direct_sam_sort_are_equivalent(self):
        k8 = ROOT / "CHARM_scripts" / "k8"
        hickit_js = ROOT / "CHARM_scripts" / "hickit.js"
        hickit = ROOT / "CHARM_scripts" / "hickit"
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            bam = work / "input.sort.bam"
            run(["samtools", "sort", "-o", str(bam), str(FIXTURE / "input.sam")])
            run(["samtools", "index", str(bam)])

            baseline = work / "baseline"
            candidate = work / "candidate"
            baseline.mkdir()
            candidate.mkdir()
            for output_dir, name_sort in (
                (
                    baseline,
                    "samtools view -h {bam} | samtools sort -n | samtools view -h".format(
                        bam=bam
                    ),
                ),
                (
                    candidate,
                    "samtools sort -n --no-PG -O SAM -T {tmp} {bam}".format(
                        tmp=work / "namesort", bam=bam
                    ),
                ),
            ):
                pipeline = (
                    "set -euo pipefail; {name_sort} | {k8} {js} sam2seg -v {snp} - "
                    "| {k8} {js} chronly - | {k8} {js} bedflt {par} - "
                    "| sed 's/-/+/g' | gzip > {output}"
                ).format(
                    name_sort=name_sort,
                    k8=k8,
                    js=hickit_js,
                    snp=FIXTURE / "phased_snps.tsv",
                    par=FIXTURE / "par.bed",
                    output=output_dir / "contacts.seg.gz",
                )
                run(["bash", "-c", pipeline])
                for distance, name in ((0, "raw.pairs.gz"), (500, "contacts.pairs.gz")):
                    with open(str(output_dir / name), "wb") as output:
                        hickit_process = subprocess.Popen(
                            [
                                str(hickit),
                                "--dup-dist=" + str(distance),
                                "-i",
                                str(output_dir / "contacts.seg.gz"),
                                "-o",
                                "-",
                            ],
                            stdout=subprocess.PIPE,
                        )
                        gzip_process = subprocess.Popen(
                            ["gzip"], stdin=hickit_process.stdout, stdout=output
                        )
                        hickit_process.stdout.close()
                        self.assertEqual(gzip_process.wait(), 0)
                        self.assertEqual(hickit_process.wait(), 0)

            old_x = int(
                run(
                    ["samtools", "view", "-c", "-q", "30", str(bam), "chrX"],
                    stdout=subprocess.PIPE,
                ).stdout
            )
            old_y = int(
                run(
                    ["samtools", "view", "-c", "-q", "30", str(bam), "chrY"],
                    stdout=subprocess.PIPE,
                ).stdout
            )
            combined = run(
                ["samtools", "view", "-q", "30", str(bam), "chrX", "chrY"],
                stdout=subprocess.PIPE,
            ).stdout.splitlines()
            new_x = sum(line.split("\t")[2] == "chrX" for line in combined)
            new_y = sum(line.split("\t")[2] == "chrY" for line in combined)
            self.assertEqual((old_x, old_y), (new_x, new_y))
            ratio = "{:.6f}\n".format(old_y / old_x)
            (baseline / "yperx.txt").write_text(ratio)
            (candidate / "yperx.txt").write_text(ratio)

            charmtools_parent = os.environ.get("CHARMTOOLS_PYTHONPATH")
            if charmtools_parent:
                environment = os.environ.copy()
                environment["PYTHONPATH"] = charmtools_parent
                for output_dir in (baseline, candidate):
                    run(
                        [
                            "python",
                            "-m",
                            "CHARMtools",
                            "clean_leg",
                            "-t",
                            "1",
                            str(output_dir / "contacts.pairs.gz"),
                            "-o",
                            str(output_dir / "c1.pairs.gz"),
                        ],
                        env=environment,
                    )
                    run(
                        [
                            "python",
                            "-m",
                            "CHARMtools",
                            "clean_isolated",
                            "-t",
                            "1",
                            "-o",
                            str(output_dir / "c12.pairs.gz"),
                            str(output_dir / "c1.pairs.gz"),
                        ],
                        env=environment,
                    )
                for name in ("c1", "c12"):
                    expected = {
                        line.rstrip("\n")
                        for line in (FIXTURE / ("expected." + name + ".pairs"))
                        .read_text()
                        .splitlines()
                    }
                    self.assertEqual(
                        set(AUDIT._lines(baseline / (name + ".pairs.gz"))), expected
                    )
            else:
                for output_dir in (baseline, candidate):
                    for name in ("c1", "c12"):
                        with gzip.open(str(output_dir / (name + ".pairs.gz")), "wt") as handle:
                            handle.write(
                                (FIXTURE / ("expected." + name + ".pairs")).read_text()
                            )

            checks = AUDIT.compare_output_dirs(baseline, candidate)
            self.assertTrue(all(passed for _, passed in checks), checks)


if __name__ == "__main__":
    unittest.main()

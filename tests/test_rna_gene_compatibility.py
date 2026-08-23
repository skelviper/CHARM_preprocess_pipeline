#!/usr/bin/env python3

import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

import pysam


SCRIPTS = Path(__file__).resolve().parents[1] / "CHARM_scripts"
sys.path.insert(0, str(SCRIPTS))

from classify_rna_gene_compatibility import (  # noqa: E402
    classify_bams,
    classify_gene_sets,
    has_confident_alignment,
    load_gene_components,
)
from summarize_rna_output_modes import summarize  # noqa: E402


SAFE_CODE = "C0123456789abcdef"


def qname(index):
    return "q{:02d}_{}_AAAAAAAA".format(index, SAFE_CODE)


def write_gtf(path):
    path.write_text(
        "".join(
            [
                'chr1\ttest\tgene\t100\t200\t.\t+\t.\tgene_id "A"; gene_name "A";\n',
                'chr1\ttest\tgene\t180\t300\t.\t-\t.\tgene_id "B"; gene_name "B";\n',
                'chr1\ttest\tgene\t400\t500\t.\t+\t.\tgene_id "C"; gene_name "C";\n',
                'chr2\ttest\tgene\t100\t200\t.\t+\t.\tgene_id "D"; gene_name "D";\n',
            ]
        )
    )


def make_record(name, genes=None):
    record = pysam.AlignedSegment()
    record.query_name = name
    record.query_sequence = "A" * 50
    record.flag = 0
    record.reference_id = 0
    record.reference_start = 100
    record.mapping_quality = 60
    record.cigar = ((0, 50),)
    record.query_qualities = pysam.qualitystring_to_array("I" * 50)
    record.set_tag("NH", 1)
    if genes:
        record.set_tag("XS", "Assigned")
        record.set_tag("XN", len(genes))
        record.set_tag("XT", ",".join(genes))
    else:
        record.set_tag("XS", "Unassigned_NoFeatures")
    return record


def write_bam(path, rows):
    header = {
        "HD": {"VN": "1.6", "SO": "queryname"},
        "SQ": [{"SN": "chr1", "LN": 1000}, {"SN": "chr2", "LN": 1000}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
        for name, genes in rows:
            handle.write(make_record(name, genes))


def bam_qnames(path):
    with pysam.AlignmentFile(str(path), "rb") as handle:
        records = list(handle.fetch(until_eof=True))
    if not all(not record.has_tag("XS") for record in records):
        raise AssertionError("filtered BAM retains XS tags")
    if not all(not record.has_tag("XN") for record in records):
        raise AssertionError("filtered BAM retains XN tags")
    if not all(not record.has_tag("XT") for record in records):
        raise AssertionError("filtered BAM retains XT tags")
    return [record.query_name for record in records]


def write_matrix(path, cells, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene"] + cells)
        writer.writerows(rows)


class RnaGeneCompatibilityTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_gene_locus_components_ignore_strand_and_preserve_disjoint_loci(self):
        annotation = self.root / "genes.gtf"
        write_gtf(annotation)
        components = load_gene_components(str(annotation))
        self.assertEqual(components["A"], components["B"])
        self.assertTrue(components["A"].isdisjoint(components["C"]))
        self.assertEqual(
            classify_gene_sets({"A"}, {"B"}, components, True), "concordant"
        )
        self.assertEqual(
            classify_gene_sets({"A"}, {"C"}, components, True), "incompatible"
        )
        self.assertEqual(
            classify_gene_sets({"A"}, set(), components, False), "r2_uninformative"
        )
        self.assertEqual(
            classify_gene_sets({"A"}, set(), components, True), "r2_genome_only"
        )
        self.assertEqual(
            classify_gene_sets(set(), {"A"}, components, True), "r1_uninformative"
        )

    def test_confident_alignment_requires_primary_unique_mapq30(self):
        confident = make_record(qname(1), None)
        low_mapq = make_record(qname(2), None)
        low_mapq.mapping_quality = 29
        multimapped = make_record(qname(3), None)
        multimapped.set_tag("NH", 2)
        self.assertTrue(has_confident_alignment([confident], 30))
        self.assertFalse(has_confident_alignment([low_mapq], 30))
        self.assertFalse(has_confident_alignment([multimapped], 30))

    def test_streaming_classifier_emits_nested_r1_cohorts_and_conserves_qnames(self):
        annotation = self.root / "genes.gtf"
        write_gtf(annotation)
        contract = self.root / "contract.json"
        contract.write_text(
            json.dumps(
                {"samples": [{"sample_name": "Cell_one", "safe_code": SAFE_CODE}]}
            )
        )
        r1_bam = self.root / "r1.bam"
        r2_bam = self.root / "r2.bam"
        write_bam(
            r1_bam,
            [
                (qname(1), ["A"]),
                (qname(2), ["A"]),
                (qname(3), ["A"]),
                (qname(4), ["A"]),
                (qname(5), None),
                (qname(6), ["A", "C"]),
                (qname(7), ["A"]),
            ],
        )
        write_bam(
            r2_bam,
            [
                (qname(1), ["A"]),
                (qname(2), ["B"]),
                (qname(3), ["C"]),
                (qname(5), ["A"]),
                (qname(6), ["C"]),
                (qname(7), None),
            ],
        )

        compatible = self.root / "compatible.bam"
        concordant = self.root / "concordant.bam"
        summary = self.root / "summary.tsv"
        classify_bams(
            str(r1_bam),
            str(r2_bam),
            str(annotation),
            str(contract),
            str(compatible),
            str(concordant),
            str(summary),
        )

        self.assertEqual(
            bam_qnames(compatible),
            [qname(1), qname(2), qname(4), qname(5), qname(6)],
        )
        self.assertEqual(bam_qnames(concordant), [qname(1), qname(2), qname(6)])
        with summary.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        cell = rows[0]
        self.assertEqual(cell["cellname"], "Cell_one")
        self.assertEqual(int(cell["r1_star_bam_qnames"]), 7)
        self.assertEqual(int(cell["r1_uninformative"]), 1)
        self.assertEqual(int(cell["r2_uninformative"]), 1)
        self.assertEqual(int(cell["r2_genome_only"]), 1)
        self.assertEqual(int(cell["concordant"]), 3)
        self.assertEqual(int(cell["incompatible"]), 1)
        self.assertEqual(int(cell["r1_compatible_qnames"]), 5)
        self.assertEqual(int(cell["r1r2_concordant_qnames"]), 3)

    def test_output_mode_summary_uses_fractional_matrix_values(self):
        modes = ["r1_all", "r1_compatible", "r1r2_concordant"]
        for mode_index, mode in enumerate(modes):
            for feature in ("gene", "exon"):
                write_matrix(
                    self.root / mode / "counts.{}.total.format.tsv".format(feature),
                    ["Cell1", "Cell2"],
                    [
                        ["G1", 3 - mode_index, 1.5],
                        ["G2", 0, 2 - mode_index / 2.0],
                    ],
                )
        per_cell = self.root / "per_cell.tsv"
        summary = self.root / "summary.tsv"
        summarize(str(self.root), modes, str(per_cell), str(summary))
        with per_cell.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(
            rows[0],
            {
                "cellname": "Cell1",
                "rna_output_type": "r1_all",
                "feature_type": "gene",
                "umi_count": "3",
                "detected_features": "1",
            },
        )
        self.assertEqual(len(rows), 12)


if __name__ == "__main__":
    unittest.main()

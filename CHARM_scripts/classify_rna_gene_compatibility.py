#!/usr/bin/env python3
"""Classify paired RNA qnames and emit filtered R1 BAM cohorts.

Inputs are query-name-sorted, featureCounts-annotated BAMs. R1 remains the
counting molecule. R2 supplies only gene-locus compatibility evidence.
"""

import argparse
import csv
import itertools
import json
import os
import re
import sys
from collections import defaultdict
from functools import cmp_to_key

import pysam


GENE_ID = re.compile(r'(?:^|;\s*)gene_id\s+"([^"]+)"')
FEATURECOUNT_TAGS = ("XS", "XN", "XT")
CATEGORIES = (
    "r1_uninformative",
    "r2_uninformative",
    "r2_genome_only",
    "concordant",
    "incompatible",
)


def compare_qnames(left, right):
    """Match samtools 1.17 strnum_cmp, including equivalent leading zeros."""
    i = j = 0
    while i < len(left) and j < len(right):
        if "0" <= left[i] <= "9" and "0" <= right[j] <= "9":
            a, b = i, j
            while i < len(left) and "0" <= left[i] <= "9":
                i += 1
            while j < len(right) and "0" <= right[j] <= "9":
                j += 1
            x, y = int(left[a:i]), int(right[b:j])
            if x != y:
                return (x > y) - (x < y)
        else:
            if left[i] != right[j]:
                return (left[i] > right[j]) - (left[i] < right[j])
            i += 1
            j += 1
    return (i < len(left)) - (j < len(right))


natural_key = cmp_to_key(compare_qnames)


def load_gene_components(gtf_path):
    spans = {}
    with open(gtf_path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError("invalid GTF row {}:{}".format(gtf_path, line_number))
            match = GENE_ID.search(fields[8])
            if match is None:
                continue
            gene_id = match.group(1)
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            key = (chrom, gene_id)
            if key in spans:
                old_start, old_end = spans[key]
                spans[key] = (min(old_start, start), max(old_end, end))
            else:
                spans[key] = (start, end)
    if not spans:
        raise ValueError("annotation contains no gene_id spans: {}".format(gtf_path))

    by_chrom = defaultdict(list)
    for (chrom, gene_id), (start, end) in spans.items():
        by_chrom[chrom].append((start, end, gene_id))

    components = defaultdict(set)
    component_number = 0
    for chrom in sorted(by_chrom):
        component_end = None
        component_token = None
        for start, end, gene_id in sorted(by_chrom[chrom]):
            if component_end is None or start > component_end:
                component_number += 1
                component_token = "{}:{}".format(chrom, component_number)
                component_end = end
            else:
                component_end = max(component_end, end)
            components[gene_id].add(component_token)
    return dict(components)


def is_confident_primary_unique(record, min_mapq):
    return not (
        record.is_unmapped
        or record.is_secondary
        or record.is_supplementary
        or record.mapping_quality < min_mapq
        or not record.has_tag("NH")
        or record.get_tag("NH") != 1
    )


def has_confident_alignment(records, min_mapq):
    return any(is_confident_primary_unique(record, min_mapq) for record in records)


def assigned_genes(records, min_mapq):
    genes = set()
    for record in records:
        if not is_confident_primary_unique(record, min_mapq):
            continue
        if not record.has_tag("XS") or record.get_tag("XS") != "Assigned":
            continue
        if not record.has_tag("XT"):
            continue
        genes.update(gene for gene in record.get_tag("XT").split(",") if gene)
    return genes


def compatible_gene_sets(r1_genes, r2_genes, components):
    if r1_genes.intersection(r2_genes):
        return True
    r1_loci = set()
    r2_loci = set()
    for gene in r1_genes:
        r1_loci.update(components.get(gene, {"gene:" + gene}))
    for gene in r2_genes:
        r2_loci.update(components.get(gene, {"gene:" + gene}))
    return bool(r1_loci.intersection(r2_loci))


def classify_gene_sets(r1_genes, r2_genes, components, r2_confidently_mapped):
    if not r1_genes:
        return "r1_uninformative"
    if not r2_genes:
        if r2_confidently_mapped:
            return "r2_genome_only"
        return "r2_uninformative"
    if compatible_gene_sets(r1_genes, r2_genes, components):
        return "concordant"
    return "incompatible"


def grouped_records(handle):
    previous_key = None
    for key, records in itertools.groupby(
        handle.fetch(until_eof=True), lambda rec: natural_key(rec.query_name)
    ):
        if previous_key is not None and key < previous_key:
            raise ValueError("BAM is not query-name sorted near {}".format(key.obj))
        previous_key = key
        # samtools can interleave distinct QNAMEs with equal numeric keys.
        by_name = defaultdict(list)
        for record in records:
            by_name[record.query_name].append(record)
        for qname in sorted(by_name):
            yield qname, by_name[qname], (key, qname)


def safe_code_from_qname(qname, allowed_codes):
    fields = qname.rsplit("_", 2)
    if len(fields) < 3 or fields[-2] not in allowed_codes:
        raise ValueError("qname lacks an authoritative safe cell code: {}".format(qname))
    return fields[-2]


def strip_featurecount_tags(records):
    for record in records:
        for tag in FEATURECOUNT_TAGS:
            if record.has_tag(tag):
                record.set_tag(tag, None)
        yield record


def load_contract(contract_path):
    with open(contract_path, encoding="utf-8") as handle:
        contract = json.load(handle)
    samples = contract.get("samples", [])
    if not samples:
        raise ValueError("input contract contains no samples")
    code_to_name = {}
    order = []
    for sample in samples:
        code = sample["safe_code"]
        name = sample["sample_name"]
        if code in code_to_name:
            raise ValueError("duplicated safe code in input contract: {}".format(code))
        code_to_name[code] = name
        order.append(code)
    return code_to_name, order


def empty_metrics():
    metrics = {category: 0 for category in CATEGORIES}
    metrics.update(
        {
            "r1_star_bam_qnames": 0,
            "r1_gene_informative": 0,
            "r2_gene_informative": 0,
            "r2_confidently_mapped": 0,
        }
    )
    return metrics


def render_summary(summary_path, metrics_by_code, code_to_name, code_order):
    fields = [
        "cellname",
        "r1_star_bam_qnames",
        "r1_gene_informative",
        "r2_gene_informative",
        "r2_confidently_mapped",
        "r1_uninformative",
        "r2_uninformative",
        "r2_genome_only",
        "concordant",
        "incompatible",
        "r1_compatible_qnames",
        "r1r2_concordant_qnames",
        "r2_genome_only_fraction_of_r1_gene_informative",
        "incompatible_fraction_of_r1_gene_informative",
        "concordant_fraction_of_r1_gene_informative",
    ]
    all_metrics = empty_metrics()
    rows = []
    for code in code_order:
        metrics = metrics_by_code[code]
        for key in all_metrics:
            all_metrics[key] += metrics[key]
        rows.append((code_to_name[code], metrics))
    rows.append(("ALL", all_metrics))

    with open(summary_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for cellname, metrics in rows:
            denominator = metrics["r1_gene_informative"]
            row = dict(metrics)
            row["cellname"] = cellname
            row["r1_compatible_qnames"] = (
                metrics["r1_star_bam_qnames"]
                - metrics["r2_genome_only"]
                - metrics["incompatible"]
            )
            row["r1r2_concordant_qnames"] = metrics["concordant"]
            row["r2_genome_only_fraction_of_r1_gene_informative"] = (
                "{:.10f}".format(metrics["r2_genome_only"] / float(denominator))
                if denominator
                else "NA"
            )
            row["incompatible_fraction_of_r1_gene_informative"] = (
                "{:.10f}".format(metrics["incompatible"] / float(denominator))
                if denominator
                else "NA"
            )
            row["concordant_fraction_of_r1_gene_informative"] = (
                "{:.10f}".format(metrics["concordant"] / float(denominator))
                if denominator
                else "NA"
            )
            writer.writerow(row)


def classify_bams(
    r1_bam,
    r2_bam,
    annotation,
    contract,
    compatible_bam,
    concordant_bam,
    summary,
    min_mapq=30,
):
    components = load_gene_components(annotation)
    code_to_name, code_order = load_contract(contract)
    metrics_by_code = {code: empty_metrics() for code in code_order}

    with pysam.AlignmentFile(r1_bam, "rb") as r1_handle, pysam.AlignmentFile(
        r2_bam, "rb"
    ) as r2_handle, pysam.AlignmentFile(
        compatible_bam, "wb", template=r1_handle
    ) as compatible_handle, pysam.AlignmentFile(
        concordant_bam, "wb", template=r1_handle
    ) as concordant_handle:
        r2_groups = iter(grouped_records(r2_handle))
        r2_current = next(r2_groups, None)
        for qname, r1_records, r1_key in grouped_records(r1_handle):
            while r2_current is not None and r2_current[2] < r1_key:
                r2_current = next(r2_groups, None)
            r2_records = []
            if r2_current is not None and r2_current[0] == qname:
                r2_records = r2_current[1]
                r2_current = next(r2_groups, None)

            r1_genes = assigned_genes(r1_records, min_mapq)
            r2_genes = assigned_genes(r2_records, min_mapq)
            r2_confidently_mapped = has_confident_alignment(r2_records, min_mapq)
            category = classify_gene_sets(
                r1_genes,
                r2_genes,
                components,
                r2_confidently_mapped,
            )
            code = safe_code_from_qname(qname, code_to_name)
            metrics = metrics_by_code[code]
            metrics["r1_star_bam_qnames"] += 1
            metrics[category] += 1
            if r1_genes:
                metrics["r1_gene_informative"] += 1
            if r2_genes:
                metrics["r2_gene_informative"] += 1
            if r2_confidently_mapped:
                metrics["r2_confidently_mapped"] += 1

            clean_records = list(strip_featurecount_tags(r1_records))
            if category not in ("r2_genome_only", "incompatible"):
                for record in clean_records:
                    compatible_handle.write(record)
            if category == "concordant":
                for record in clean_records:
                    concordant_handle.write(record)

    for code, metrics in metrics_by_code.items():
        categorized = sum(metrics[category] for category in CATEGORIES)
        if categorized != metrics["r1_star_bam_qnames"]:
            raise ValueError("qname category conservation failed for {}".format(code))
    render_summary(summary, metrics_by_code, code_to_name, code_order)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1-bam", required=True)
    parser.add_argument("--r2-bam", required=True)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--contract", required=True)
    parser.add_argument("--compatible-bam", required=True)
    parser.add_argument("--concordant-bam", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--min-mapq", type=int, default=30)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        classify_bams(
            args.r1_bam,
            args.r2_bam,
            args.annotation,
            args.contract,
            args.compatible_bam,
            args.concordant_bam,
            args.summary,
            args.min_mapq,
        )
    except (OSError, ValueError, pysam.SamtoolsError) as error:
        print("RNA gene-compatibility classification error: {}".format(error), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

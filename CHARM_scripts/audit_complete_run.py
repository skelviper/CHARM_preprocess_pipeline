#!/usr/bin/env python3
"""Audit the declared CHARM/HiRES deliverables after metadata generation."""

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path


FEATURES = ("gene", "exon")
AUDIT_OUTPUT = "qc/COMPLETE_RUN_AUDIT.tsv"
COMPATIBILITY_COLUMNS = {
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
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_matrix(path, expected_cells):
    sums = {cell: 0.0 for cell in expected_cells}
    detected = {cell: 0 for cell in expected_cells}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header != ["gene"] + expected_cells:
            raise ValueError("matrix header mismatch: {}".format(path))
        for line_number, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise ValueError(
                    "matrix width mismatch: {}:{}".format(path, line_number)
                )
            for cell, value in zip(expected_cells, row[1:]):
                numeric = float(value)
                if not math.isfinite(numeric) or numeric < 0:
                    raise ValueError(
                        "invalid matrix value: {}:{}".format(path, line_number)
                    )
                sums[cell] += numeric
                detected[cell] += int(numeric != 0)
    return sums, detected


def require_numeric(row, columns, cell):
    for column in columns:
        value = float(row[column])
        if not math.isfinite(value) or value < 0:
            raise ValueError("invalid metadata value {} for {}".format(column, cell))


def read_pair_stat(path, expected_cells):
    suffix = ".pairs.gz"
    observed = []
    values = {}
    with path.open(newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(row) != 2:
                raise ValueError(
                    "pair statistic width mismatch: {}:{}".format(path, line_number)
                )
            basename = Path(row[0]).name
            if not basename.endswith(suffix):
                raise ValueError("invalid pair statistic path: {}".format(row[0]))
            cell = basename[: -len(suffix)]
            try:
                value = int(row[1])
            except ValueError:
                raise ValueError("invalid pair statistic count: {}".format(row[1]))
            if value < 0 or cell in values:
                raise ValueError("invalid pair statistic cell: {}".format(cell))
            observed.append(cell)
            values[cell] = value
    if observed != expected_cells:
        raise ValueError("pair statistic cell order or membership mismatch: {}".format(path))
    return values


def normalized_contract_hash(samples):
    normalized = []
    for sample in samples:
        row = {
            "sample_name": sample["sample_name"],
            "safe_code": sample["safe_code"],
        }
        for mate in ("r1", "r2"):
            read = sample["reads"][mate]
            row["{}_bytes".format(mate)] = read["target_stat"]["size"]
            row["{}_first_record_sha256".format(mate)] = read[
                "first_record_sha256"
            ]
        normalized.append(row)
    payload = json.dumps(normalized, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def audit(work_dir, pipeline_dir):
    checks = []

    contract_path = work_dir / "qc/input_contract/current.json"
    with contract_path.open() as handle:
        contract = json.load(handle)
    samples = contract["samples"]
    cells = [sample["sample_name"] for sample in samples]
    if not cells or len(cells) != contract["sample_count"]:
        raise ValueError("input contract sample count mismatch")
    if len(set(cells)) != len(cells):
        raise ValueError("input contract contains duplicate cells")
    checks.append(("input_contract", "PASS", "{} unique cells".format(len(cells))))
    checks.append(
        ("normalized_contract_sha256", "PASS", normalized_contract_hash(samples))
    )

    config_path = work_dir / "qc/provenance/effective_config.json"
    with config_path.open() as handle:
        config = json.load(handle)
    modes = tuple(config["rna_output_types"])
    if not modes or len(set(modes)) != len(modes):
        raise ValueError("invalid effective RNA output modes")
    experiment_type = config["experiment_type"]
    if experiment_type not in ("charm", "hires"):
        raise ValueError("invalid effective experiment type")
    primary_mode = config["rna_primary_output_type"]
    if primary_mode not in modes:
        raise ValueError("primary RNA mode is not selected")

    target_path = work_dir / "qc/target_outputs.tsv"
    with target_path.open(newline="") as handle:
        targets = list(csv.DictReader(handle, delimiter="\t"))
    all_outputs = [
        work_dir / row["output"]
        for row in targets
        if row["target"] == "all"
        and row["enabled"] == "1"
        and row["output"] not in ("NA", AUDIT_OUTPUT)
    ]
    missing = [
        str(path.relative_to(work_dir)) for path in all_outputs if not path.is_file()
    ]
    empty = [
        str(path.relative_to(work_dir))
        for path in all_outputs
        if path.is_file() and path.stat().st_size == 0
    ]
    if missing or empty:
        raise ValueError(
            "all-target failure; missing={!r}; empty={!r}".format(missing, empty)
        )
    checks.append(
        ("all_target_outputs", "PASS", "{} nonempty files".format(len(all_outputs)))
    )

    matrix_metrics = {}
    for mode in modes:
        for feature in FEATURES:
            path = work_dir / "result/RNA_Res/{}/counts.{}.total.format.tsv".format(
                mode, feature
            )
            sums, detected = read_matrix(path, cells)
            matrix_metrics[(mode, feature)] = {
                "umis": sums,
                "detected": detected,
            }
    checks.append(
        (
            "rna_matrix_contract",
            "PASS",
            "{} authoritative matrices".format(len(modes) * len(FEATURES)),
        )
    )

    mode_table_path = work_dir / "qc/stat/rna.output_modes.per_cell.tsv"
    with mode_table_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != [
            "cellname",
            "rna_output_type",
            "feature_type",
            "umi_count",
            "detected_features",
        ]:
            raise ValueError("RNA mode per-cell table header mismatch")
        mode_rows = list(reader)
    expected_mode_keys = {
        (cell, mode, feature)
        for mode in modes
        for feature in FEATURES
        for cell in cells
    }
    observed_mode_keys = set()
    for row in mode_rows:
        key = (row["cellname"], row["rna_output_type"], row["feature_type"])
        if key in observed_mode_keys:
            raise ValueError("duplicate RNA mode per-cell row: {}".format(key))
        observed_mode_keys.add(key)
        if key not in expected_mode_keys:
            raise ValueError("unexpected RNA mode per-cell row: {}".format(key))
        cell, mode, feature = key
        metrics = matrix_metrics[(mode, feature)]
        if abs(float(row["umi_count"]) - metrics["umis"][cell]) > 1e-6:
            raise ValueError("RNA mode per-cell UMI mismatch: {}".format(key))
        if int(row["detected_features"]) != metrics["detected"][cell]:
            raise ValueError(
                "RNA mode per-cell detected-feature mismatch: {}".format(key)
            )
    if observed_mode_keys != expected_mode_keys:
        raise ValueError("RNA mode per-cell table is incomplete")

    summary_path = work_dir / "qc/stat/rna.output_modes.summary.tsv"
    with summary_path.open(newline="") as handle:
        summary_rows = list(csv.DictReader(handle, delimiter="\t"))
    expected_summary_keys = {(mode, feature) for mode in modes for feature in FEATURES}
    observed_summary_keys = set()
    for row in summary_rows:
        key = (row["rna_output_type"], row["feature_type"])
        if key in observed_summary_keys or key not in expected_summary_keys:
            raise ValueError("invalid RNA mode summary row: {}".format(key))
        observed_summary_keys.add(key)
        metrics = matrix_metrics[key]
        per_cell_umis = [metrics["umis"][cell] for cell in cells]
        expected_values = {
            "cells": len(cells),
            "total_umis": sum(per_cell_umis),
            "median_umis_per_cell": statistics.median(per_cell_umis),
            "min_umis_per_cell": min(per_cell_umis),
            "max_umis_per_cell": max(per_cell_umis),
            "total_detected_cell_features": sum(metrics["detected"].values()),
        }
        for column, expected in expected_values.items():
            if abs(float(row[column]) - expected) > 1e-6:
                raise ValueError(
                    "RNA mode summary mismatch: {} {}".format(key, column)
                )
    if observed_summary_keys != expected_summary_keys:
        raise ValueError("RNA mode summary table is incomplete")
    checks.append(
        (
            "rna_mode_tables",
            "PASS",
            "{} per-cell rows; {} summary rows".format(
                len(mode_rows), len(summary_rows)
            ),
        )
    )

    if any(mode != "r1_all" for mode in modes):
        compatibility_path = work_dir / "qc/stat/rna.gene_compatibility.per_cell.tsv"
        with compatibility_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if [row["cellname"] for row in rows] != cells + ["ALL"]:
            raise ValueError("RNA compatibility row contract mismatch")
        categories = (
            "r1_uninformative",
            "r2_uninformative",
            "r2_genome_only",
            "concordant",
            "incompatible",
        )
        for row in rows:
            total = int(row["r1_star_bam_qnames"])
            partition = sum(int(row[column]) for column in categories)
            compatible = (
                total - int(row["r2_genome_only"]) - int(row["incompatible"])
            )
            if partition != total or compatible != int(row["r1_compatible_qnames"]):
                raise ValueError(
                    "RNA compatibility conservation failed: {}".format(
                        row["cellname"]
                    )
                )
            if int(row["concordant"]) != int(row["r1r2_concordant_qnames"]):
                raise ValueError(
                    "RNA concordant nesting failed: {}".format(row["cellname"])
                )
        checks.append(("rna_qname_conservation", "PASS", "cells plus ALL"))

    metadata_path = work_dir / "qc/metadata_raw.tsv"
    with metadata_path.open(newline="") as handle:
        metadata = list(csv.DictReader(handle, delimiter="\t"))
    if [row["cellname"] for row in metadata] != cells:
        raise ValueError("metadata cell order or membership mismatch")
    required_columns = {
        "cellname",
        "experiment_type",
        "RNA_primary_output_type",
        "UMIs_gene",
        "genes_gene",
        "UMIs_exon",
        "genes_exon",
        "raw_pairs",
        "pairs_dedup",
        "pairs_clean1",
        "pairs_clean2",
        "rna_total_mapped_reads",
        "rna_assigned_reads",
        "rna_annotation_rate",
        "rna_dedup_rate",
        "rna_clean_reads",
        "rna_gatc_reads",
        "rna_dna_contam_rate",
        "pairs_clean3",
        "inter_pairs_clean3",
        "pairsValidRatio",
        "interPairsRatio",
    }
    if experiment_type == "charm":
        required_columns.update(
            {
                "atac_reads",
                "ct_reads",
                "atac_dedup_rate",
                "ct_dedup_rate",
                "nCount_atac",
                "nCount_ct",
                "TSS.enrichment.atac",
                "TSS.enrichment.ct",
            }
        )
    observed_columns = set(metadata[0]) if metadata else set()
    if not required_columns.issubset(observed_columns):
        raise ValueError(
            "missing metadata columns: {}".format(
                sorted(required_columns - observed_columns)
            )
        )
    mode_specific_columns = {
        "{}_{}_{}".format(prefix, feature, mode)
        for prefix in ("UMIs", "genes")
        for feature in FEATURES
        for mode in modes
    }
    non_primary_columns = mode_specific_columns | COMPATIBILITY_COLUMNS | {
        "UMIs_gene_genome1",
        "genes_gene_genome1",
        "UMIs_gene_genome2",
        "genes_gene_genome2",
    }
    unexpected_columns = sorted(observed_columns & non_primary_columns)
    if unexpected_columns:
        raise ValueError(
            "non-primary RNA columns leaked into metadata: {}".format(
                unexpected_columns
            )
        )
    text_columns = {"cellname", "experiment_type", "RNA_primary_output_type"}
    numeric_columns = sorted(required_columns - text_columns)
    c123_counts = read_pair_stat(work_dir / "qc/stat/pairs.c123.stat", cells)
    c123_inter_counts = read_pair_stat(
        work_dir / "qc/stat/inter.pairs.c123.stat", cells
    )
    for row in metadata:
        if (
            row["experiment_type"] != experiment_type
            or row["RNA_primary_output_type"] != primary_mode
        ):
            raise ValueError("metadata experiment or primary RNA mode mismatch")
        require_numeric(row, numeric_columns, row["cellname"])
        cell = row["cellname"]
        if int(float(row["pairs_clean3"])) != c123_counts[cell]:
            raise ValueError("metadata/c123 count mismatch: {}".format(cell))
        if int(float(row["inter_pairs_clean3"])) != c123_inter_counts[cell]:
            raise ValueError("metadata/c123 inter-count mismatch: {}".format(cell))
        raw_pairs = float(row["raw_pairs"])
        expected_valid_ratio = c123_counts[cell] / raw_pairs if raw_pairs > 0 else 0
        expected_inter_ratio = (
            c123_inter_counts[cell] / c123_counts[cell]
            if c123_counts[cell] > 0
            else 0
        )
        if abs(float(row["pairsValidRatio"]) - expected_valid_ratio) > 1e-8:
            raise ValueError("metadata/c123 valid-ratio mismatch: {}".format(cell))
        if abs(float(row["interPairsRatio"]) - expected_inter_ratio) > 1e-8:
            raise ValueError("metadata/c123 inter-ratio mismatch: {}".format(cell))
        for feature in FEATURES:
            umi_column = "UMIs_{}".format(feature)
            detected_column = "genes_{}".format(feature)
            metrics = matrix_metrics[(primary_mode, feature)]
            if abs(float(row[umi_column]) - metrics["umis"][row["cellname"]]) > 1e-6:
                raise ValueError(
                    "metadata/matrix mismatch: {} {}".format(
                        row["cellname"], umi_column
                    )
                )
            if int(float(row[detected_column])) != metrics["detected"][row["cellname"]]:
                raise ValueError(
                    "metadata/matrix mismatch: {} {}".format(
                        row["cellname"], detected_column
                    )
                )
    checks.append(
        (
            "metadata_contract",
            "PASS",
            "{} rows; {} columns".format(len(metadata), len(observed_columns)),
        )
    )

    provenance_path = work_dir / "qc/provenance/source_files.sha256.tsv"
    with provenance_path.open(newline="") as handle:
        provenance = list(csv.DictReader(handle, delimiter="\t"))
    mismatches = []
    for row in provenance:
        source = (pipeline_dir / row["path"]).resolve()
        if not source.is_file() or sha256(source) != row["sha256"]:
            mismatches.append(row["path"])
    if mismatches:
        raise ValueError("source provenance mismatch: {}".format(mismatches))
    checks.append(
        ("source_provenance", "PASS", "{} hashes".format(len(provenance)))
    )

    receipt_path = work_dir / AUDIT_OUTPUT
    with receipt_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("check", "status", "detail"))
        writer.writerows(checks)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--pipeline-dir", required=True, type=Path)
    args = parser.parse_args()
    audit(args.work_dir.resolve(), args.pipeline_dir.resolve())


if __name__ == "__main__":
    main()

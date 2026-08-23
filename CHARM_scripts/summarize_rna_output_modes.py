#!/usr/bin/env python3
"""Summarize mode-specific RNA matrices per cell and across the cohort."""

import argparse
import csv
import os
import statistics
import sys


FEATURES = ("gene", "exon")


def read_matrix(path):
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("empty RNA matrix: {}".format(path))
        if not header or header[0] != "gene":
            raise ValueError("RNA matrix first column must be gene: {}".format(path))
        cells = header[1:]
        totals = [0.0] * len(cells)
        detected = [0] * len(cells)
        for row_number, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise ValueError("matrix row {} has wrong width: {}".format(row_number, path))
            for index, value in enumerate(row[1:]):
                number = float(value)
                totals[index] += number
                detected[index] += int(number != 0)
    return cells, totals, detected


def format_number(value):
    return "{:.10f}".format(value).rstrip("0").rstrip(".") or "0"


def summarize(matrix_root, modes, per_cell_path, summary_path):
    per_cell_rows = []
    summary_rows = []
    expected_cells = None
    for mode in modes:
        for feature in FEATURES:
            path = os.path.join(matrix_root, mode, "counts.{}.total.format.tsv".format(feature))
            cells, totals, detected = read_matrix(path)
            if expected_cells is None:
                expected_cells = cells
            elif cells != expected_cells:
                raise ValueError("cell columns differ between RNA output modes")
            for cell, total, detected_features in zip(cells, totals, detected):
                per_cell_rows.append(
                    {
                        "cellname": cell,
                        "rna_output_type": mode,
                        "feature_type": feature,
                        "umi_count": format_number(total),
                        "detected_features": detected_features,
                    }
                )
            summary_rows.append(
                {
                    "rna_output_type": mode,
                    "feature_type": feature,
                    "cells": len(cells),
                    "total_umis": format_number(sum(totals)),
                    "median_umis_per_cell": format_number(statistics.median(totals)),
                    "min_umis_per_cell": format_number(min(totals) if totals else 0),
                    "max_umis_per_cell": format_number(max(totals) if totals else 0),
                    "total_detected_cell_features": sum(detected),
                }
            )

    os.makedirs(os.path.dirname(os.path.abspath(per_cell_path)), exist_ok=True)
    with open(per_cell_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(per_cell_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(per_cell_rows)
    with open(summary_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(summary_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix-root", required=True)
    parser.add_argument("--modes", nargs="+", required=True)
    parser.add_argument("--per-cell", required=True)
    parser.add_argument("--summary", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        summarize(args.matrix_root, args.modes, args.per_cell, args.summary)
    except (OSError, ValueError) as error:
        print("RNA output-mode summary error: {}".format(error), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

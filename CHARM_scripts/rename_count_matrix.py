#!/usr/bin/env python3
"""Restore authoritative sample names in an UMI-tools count matrix.

This keeps the safe-code/name mapping approach validated by Part 2
``rna_read_mode/scripts/rename_count_matrix.py``. UMI-tools omits cells with no
counted UMIs, so expected-but-missing columns are restored as zeroes. Unknown or
duplicated cell codes remain contract errors.
"""

import argparse
import csv
import os
import sys
import tempfile

from input_contract import InputContractError, load_contract_file, write_if_changed


def rename_count_matrix(input_path, contract_path, output_path, receipt_path):
    contract = load_contract_file(contract_path)
    samples = contract["samples"]
    expected_codes = [sample["safe_code"] for sample in samples]
    code_to_name = {
        sample["safe_code"]: sample["sample_name"] for sample in samples
    }

    try:
        source_handle = open(input_path, "r", newline="")
    except OSError as error:
        raise ValueError("cannot open UMI count matrix: {}".format(error))

    output_dir = os.path.dirname(os.path.abspath(output_path))
    os.makedirs(output_dir, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=".rename-matrix.", dir=output_dir)
    row_count = 0
    try:
        with source_handle, os.fdopen(descriptor, "w", newline="") as destination:
            reader = csv.reader(source_handle, delimiter="\t")
            writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
            try:
                header = next(reader)
            except StopIteration:
                raise ValueError("empty UMI count matrix: {}".format(input_path))
            if not header or header[0] != "gene":
                raise ValueError("UMI count matrix first column must be 'gene'")
            observed_codes = header[1:]
            if len(observed_codes) != len(set(observed_codes)):
                duplicates = sorted(
                    code for code in set(observed_codes) if observed_codes.count(code) > 1
                )
                raise ValueError(
                    "duplicated cell codes in UMI count matrix: {}".format(
                        ",".join(duplicates)
                    )
                )
            observed_code_set = set(observed_codes)
            missing = [
                code for code in expected_codes if code not in observed_code_set
            ]
            unexpected = sorted(observed_code_set - set(expected_codes))
            if unexpected:
                raise ValueError(
                    "unexpected cell codes in UMI count matrix: {}".format(
                        ",".join(unexpected)
                    )
                )

            positions = {
                code: position + 1 for position, code in enumerate(observed_codes)
            }
            writer.writerow(["gene"] + [code_to_name[code] for code in expected_codes])
            for row_number, row in enumerate(reader, start=2):
                if len(row) != len(header):
                    raise ValueError(
                        "UMI count matrix row {} has {} columns, expected {}".format(
                            row_number, len(row), len(header)
                        )
                    )
                writer.writerow(
                    [row[0]]
                    + [
                        row[positions[code]] if code in positions else "0"
                        for code in expected_codes
                    ]
                )
                row_count += 1
            destination.flush()
            os.fsync(destination.fileno())
        os.replace(temporary, output_path)
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise

    receipt = (
        "metric\tvalue\n"
        "contract_sha256\t{}\n"
        "expected_cells\t{}\n"
        "observed_cells\t{}\n"
        "missing_cells_filled_zero\t{}\n"
        "unexpected_cells\t{}\n"
        "missing_cell_codes\t{}\n"
        "missing_cell_names\t{}\n"
        "feature_rows\t{}\n"
        "cell_column_contract\tPASS\n"
    ).format(
        contract["contract_sha256"],
        len(expected_codes),
        len(observed_codes),
        len(missing),
        len(unexpected),
        ",".join(missing) or "none",
        ",".join(code_to_name[code] for code in missing) or "none",
        row_count,
    )
    write_if_changed(receipt_path, receipt.encode("utf-8"))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--contract", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--receipt", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        rename_count_matrix(
            args.input, args.contract, args.output, args.receipt
        )
    except (InputContractError, OSError, ValueError) as error:
        print("count-matrix rename error: {}".format(error), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

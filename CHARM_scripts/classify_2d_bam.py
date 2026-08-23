#!/usr/bin/env python3
"""Classify a DNA BAM before the CHARM BAM-to-SEG conversion.

The machine-readable TSV is authoritative.  The legacy yperx scalar is kept
only for consumers that still require a finite decimal; the receipt records
whether that scalar is measured or operational.
"""

import argparse
import csv
from decimal import Decimal, InvalidOperation
import os
from pathlib import Path
import subprocess
import sys
import tempfile


CONTRACT_VERSION = "2d_bam_sex_v1"
FIELDS = (
    "contract_version",
    "input_state",
    "sex_state",
    "sex_call",
    "no_xy_fallback",
    "fallback_applied",
    "ratio_state",
    "measured",
    "numerator_y",
    "denominator_x",
    "measured_yperx",
    "effective_yperx",
    "total_records",
    "mapped_records",
    "qualifying_records",
    "x_records",
    "y_records",
    "contact_min_mapq",
    "sex_min_mapq",
    "yperx_threshold",
)


class ClassificationError(ValueError):
    pass


def _finite_decimal(value, label):
    try:
        result = Decimal(str(value))
    except InvalidOperation:
        raise ClassificationError("{} is not a decimal: {!r}".format(label, value))
    if not result.is_finite():
        raise ClassificationError("{} must be finite: {!r}".format(label, value))
    return result


def _format_ratio(numerator, denominator):
    return "{:.6f}".format(float(numerator) / float(denominator))


def scan_alignment_records(bam_path, samtools, contact_min_mapq, sex_min_mapq):
    counts = {
        "total_records": 0,
        "mapped_records": 0,
        "qualifying_records": 0,
        "x_records": 0,
        "y_records": 0,
    }
    with tempfile.TemporaryFile(mode="w+t") as stderr_handle:
        try:
            process = subprocess.Popen(
                [samtools, "view", str(bam_path)],
                stdout=subprocess.PIPE,
                stderr=stderr_handle,
                universal_newlines=True,
            )
        except OSError as error:
            raise ClassificationError("cannot execute samtools: {}".format(error))

        try:
            for line_number, line in enumerate(process.stdout, 1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 11:
                    raise ClassificationError(
                        "samtools emitted malformed SAM at record {}".format(line_number)
                    )
                try:
                    flag = int(fields[1])
                    mapq = int(fields[4])
                except ValueError:
                    raise ClassificationError(
                        "samtools emitted invalid FLAG/MAPQ at record {}".format(
                            line_number
                        )
                    )
                counts["total_records"] += 1
                if flag & 0x4:
                    continue
                counts["mapped_records"] += 1
                if mapq >= contact_min_mapq:
                    counts["qualifying_records"] += 1
                if mapq >= sex_min_mapq:
                    if fields[2] == "chrX":
                        counts["x_records"] += 1
                    elif fields[2] == "chrY":
                        counts["y_records"] += 1
        except Exception:
            process.terminate()
            process.wait()
            raise

        return_code = process.wait()
        stderr_handle.seek(0)
        stderr_text = stderr_handle.read().strip()
    if return_code != 0:
        detail = ": {}".format(stderr_text) if stderr_text else ""
        raise ClassificationError(
            "samtools view failed with exit code {}{}".format(return_code, detail)
        )
    return counts


def classify_counts(
    counts,
    contact_min_mapq,
    sex_min_mapq,
    yperx_threshold,
    no_xy_fallback,
):
    threshold = _finite_decimal(yperx_threshold, "yperx threshold")
    if threshold < 0:
        raise ClassificationError("yperx threshold must be nonnegative")
    if contact_min_mapq < 0 or sex_min_mapq < 0:
        raise ClassificationError("MAPQ thresholds must be nonnegative")
    fallback = no_xy_fallback.upper()
    if fallback not in ("XX", "XY"):
        raise ClassificationError("no-XY fallback must be XX or XY")

    row = {
        "contract_version": CONTRACT_VERSION,
        "input_state": "qualifying_mapped_records",
        "sex_state": "not_applicable",
        "sex_call": "",
        "no_xy_fallback": fallback,
        "fallback_applied": "false",
        "ratio_state": "not_applicable",
        "measured": "false",
        "numerator_y": str(counts["y_records"]),
        "denominator_x": str(counts["x_records"]),
        "measured_yperx": "",
        "effective_yperx": "",
        "contact_min_mapq": str(contact_min_mapq),
        "sex_min_mapq": str(sex_min_mapq),
        "yperx_threshold": str(yperx_threshold),
    }
    row.update({key: str(value) for key, value in counts.items()})

    if counts["total_records"] == 0:
        row["input_state"] = "header_only"
        return row
    if counts["mapped_records"] == 0:
        row["input_state"] = "unmapped_only"
        return row
    if counts["qualifying_records"] == 0:
        row["input_state"] = "low_mapq_only"
        return row

    x_records = counts["x_records"]
    y_records = counts["y_records"]
    if x_records == 0 and y_records == 0:
        row["sex_state"] = "undetermined_no_xy"
        row["sex_call"] = fallback
        row["fallback_applied"] = "true"
        row["ratio_state"] = "undefined_zero_zero"
        if fallback == "XX":
            row["effective_yperx"] = "0.000000"
        else:
            row["effective_yperx"] = "{:.6f}".format(float(threshold + 1))
        return row

    if x_records == 0:
        row["sex_state"] = "xy_y_only"
        row["sex_call"] = "XY"
        row["ratio_state"] = "positive_infinity"
        # A denominator guard keeps the legacy scalar finite.  The receipt
        # remains authoritative and marks this value as not measured.
        row["effective_yperx"] = _format_ratio(y_records, 1)
        return row

    measured_yperx = _format_ratio(y_records, x_records)
    sex_call = "XY" if Decimal(measured_yperx) > threshold else "XX"
    row["measured"] = "true"
    row["ratio_state"] = "finite"
    row["measured_yperx"] = measured_yperx
    row["effective_yperx"] = measured_yperx
    row["sex_call"] = sex_call
    if y_records == 0:
        row["sex_state"] = "xx_x_only"
    else:
        row["sex_state"] = "{}_ratio".format(sex_call.lower())
    return row


def _atomic_write(path, payload):
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=destination.name + ".", suffix=".tmp", dir=str(destination.parent)
    )
    try:
        with os.fdopen(descriptor, "w") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, str(destination))
    finally:
        if os.path.exists(temporary_name):
            os.unlink(temporary_name)


def write_outputs(row, classification_output, yperx_output):
    from io import StringIO

    buffer = StringIO()
    writer = csv.DictWriter(
        buffer, fieldnames=FIELDS, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerow(row)
    _atomic_write(classification_output, buffer.getvalue())
    legacy_scalar = row["effective_yperx"]
    _atomic_write(yperx_output, legacy_scalar + "\n" if legacy_scalar else "")


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--classification-output", required=True)
    parser.add_argument("--yperx-output", required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--contact-min-mapq", type=int, default=20)
    parser.add_argument("--sex-min-mapq", type=int, default=30)
    parser.add_argument("--yperx-threshold", required=True)
    parser.add_argument("--no-xy-fallback", choices=("XX", "XY"), default="XX")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        counts = scan_alignment_records(
            args.bam, args.samtools, args.contact_min_mapq, args.sex_min_mapq
        )
        row = classify_counts(
            counts,
            args.contact_min_mapq,
            args.sex_min_mapq,
            args.yperx_threshold,
            args.no_xy_fallback,
        )
        write_outputs(row, args.classification_output, args.yperx_output)
    except ClassificationError as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        return 1
    print(
        "|".join(
            (
                row["input_state"],
                row["sex_state"],
                row["sex_call"] or "NA",
                row["qualifying_records"],
            )
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

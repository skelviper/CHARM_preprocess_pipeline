#!/usr/bin/env python3
"""Contract-bound helpers for legacy CHARM statistics generation.

The worker reads the launcher-frozen input contract without rediscovering
``Rawdata``.  It provides the exact cohort and safe-code mapping to the legacy
statistics wrapper and fails closed on malformed RNA read identifiers or
missing enabled structure outputs.
"""

import argparse
import csv
import gzip
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from input_contract import InputContractError, load_contract_file, write_if_changed


SAFE_CODE_RE = re.compile(r"^C[0-9a-f]{16}$")
UMI_RE = re.compile(r"^[ACGTNacgtn]{8}$")
RMSD_RE = re.compile(
    r"top3 RMS RMSD:\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
)
MANIFEST_FIELDS = ["sample", "safe_code"]


class StatisticsContractError(ValueError):
    pass


def load_frozen_contract(contract_path, work_dir):
    """Load the signed contract only; never enumerate the live Rawdata tree."""

    contract = load_contract_file(contract_path)
    expected_work_dir = os.path.abspath(work_dir)
    if os.path.abspath(contract.get("work_dir", "")) != expected_work_dir:
        raise StatisticsContractError(
            "input contract belongs to {!r}, not {!r}".format(
                contract.get("work_dir"), expected_work_dir
            )
        )

    seen_codes = {}
    for sample in contract["samples"]:
        sample_name = sample.get("sample_name")
        safe_code = sample.get("safe_code")
        if not isinstance(sample_name, str) or not sample_name:
            raise StatisticsContractError("input contract has an invalid sample name")
        if not isinstance(safe_code, str) or not SAFE_CODE_RE.fullmatch(safe_code):
            raise StatisticsContractError(
                "input contract has malformed safe code for {!r}: {!r}".format(
                    sample_name, safe_code
                )
            )
        previous = seen_codes.get(safe_code)
        if previous is not None:
            raise StatisticsContractError(
                "input contract maps safe code {} to both {!r} and {!r}".format(
                    safe_code, previous, sample_name
                )
            )
        seen_codes[safe_code] = sample_name

        reads = sample.get("reads")
        if not isinstance(reads, dict) or set(reads) != {"r1", "r2"}:
            raise StatisticsContractError(
                "input contract has invalid read entries for {!r}".format(sample_name)
            )
        for mate in ("r1", "r2"):
            logical = reads[mate].get("logical_path")
            resolved = reads[mate].get("resolved_path")
            if not isinstance(logical, str) or os.path.isabs(logical):
                raise StatisticsContractError(
                    "input contract has invalid {} logical path for {!r}".format(
                        mate, sample_name
                    )
                )
            if not isinstance(resolved, str) or not os.path.isabs(resolved):
                raise StatisticsContractError(
                    "input contract has invalid {} resolved path for {!r}".format(
                        mate, sample_name
                    )
                )
    return contract


def manifest_rows(contract):
    rows = []
    for sample in contract["samples"]:
        rows.append(
            {
                "sample": sample["sample_name"],
                "safe_code": sample["safe_code"],
            }
        )
    return rows


def _tsv_bytes(rows, fields):
    from io import StringIO

    buffer = StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def write_manifest(contract, output_path):
    write_if_changed(output_path, _tsv_bytes(manifest_rows(contract), MANIFEST_FIELDS))


def _gzip_newline_count(path):
    count = 0
    try:
        with gzip.open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                count += chunk.count(b"\n")
    except (EOFError, OSError) as error:
        raise StatisticsContractError(
            "cannot count frozen FASTQ {!r}: {}".format(path, error)
        )
    return count


def write_raw_fastq_stats(contract, output_path, workers):
    samples = contract["samples"]
    resolved_paths = [sample["reads"]["r1"]["resolved_path"] for sample in samples]
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        counts = list(pool.map(_gzip_newline_count, resolved_paths))
    lines = []
    for sample, count in zip(samples, counts):
        logical = sample["reads"]["r1"]["logical_path"]
        label = logical if logical.startswith("./") else "./" + logical
        lines.append("{}\t{}\n".format(label, count))
    write_if_changed(output_path, "".join(lines).encode("utf-8"))


def code_to_sample(contract):
    mapping = {}
    for sample in contract["samples"]:
        code = sample["safe_code"]
        if code in mapping:
            raise StatisticsContractError(
                "duplicate safe code in input contract: {}".format(code)
            )
        mapping[code] = sample["sample_name"]
    return mapping


def sample_from_qname(qname, mapping):
    try:
        prefix, code, umi = qname.rsplit("_", 2)
    except ValueError:
        raise StatisticsContractError(
            "RNA QNAME lacks the documented penultimate safe-code token: {!r}".format(
                qname
            )
        )
    if not prefix or not SAFE_CODE_RE.fullmatch(code) or not UMI_RE.fullmatch(umi):
        raise StatisticsContractError(
            "malformed RNA QNAME cell/UMI suffix: {!r}".format(qname)
        )
    try:
        return mapping[code]
    except KeyError:
        raise StatisticsContractError(
            "RNA QNAME contains unknown safe code {}".format(code)
        )


def write_rna_alignment_stats(contract, input_handle, output_path):
    mapping = code_to_sample(contract)
    total = {}
    assigned = {}
    for line_number, line in enumerate(input_handle, start=1):
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 11:
            raise StatisticsContractError(
                "malformed SAM row {}: expected at least 11 fields".format(line_number)
            )
        sample = sample_from_qname(fields[0], mapping)
        total[sample] = total.get(sample, 0) + 1
        if any(field.startswith("XS:Z:Assigned") for field in fields[11:]):
            assigned[sample] = assigned.get(sample, 0) + 1
    lines = [
        "{}\t{}\t{}\n".format(sample, total[sample], assigned.get(sample, 0))
        for sample in sorted(total)
    ]
    write_if_changed(output_path, "".join(lines).encode("utf-8"))


def _gzip_pair_counts(path):
    rows = 0
    inter = 0
    try:
        with gzip.open(path, "rt") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 4:
                    raise StatisticsContractError(
                        "malformed pairs row {} in {!r}".format(line_number, path)
                    )
                rows += 1
                inter += fields[1] != fields[3]
    except (EOFError, OSError) as error:
        raise StatisticsContractError(
            "cannot read required c123 pairs {!r}: {}".format(path, error)
        )
    return rows, inter


def write_cleaned_pair_stats(contract, work_dir, pairs_output, inter_output, workers):
    work_path = Path(work_dir)
    samples = [sample["sample_name"] for sample in contract["samples"]]
    pair_paths = [
        work_path / "result" / "cleaned_pairs" / "c123" / "{}.pairs.gz".format(sample)
        for sample in samples
    ]
    missing_pairs = [str(path) for path in pair_paths if not path.is_file()]
    if missing_pairs:
        raise StatisticsContractError(
            "2D c123 QC is missing required output(s): {}".format(
                ", ".join(missing_pairs)
            )
        )

    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        counts = list(pool.map(_gzip_pair_counts, pair_paths))

    pair_lines = []
    inter_lines = []
    for path, (rows, inter) in zip(pair_paths, counts):
        label = "./" + str(path.relative_to(work_path))
        pair_lines.append("{}\t{}\n".format(label, rows))
        inter_lines.append("{}\t{}\n".format(label, inter))

    write_if_changed(pairs_output, "".join(pair_lines).encode("utf-8"))
    write_if_changed(inter_output, "".join(inter_lines).encode("utf-8"))


def write_rmsd_stats(contract, work_dir, resolutions, output_path):
    work_path = Path(work_dir)
    samples = [sample["sample_name"] for sample in contract["samples"]]
    rmsd_paths = [
        (
            sample,
            resolution,
            work_path
            / "result"
            / "3d_info"
            / sample
            / "{}.{}.align.rms.info".format(sample, resolution),
        )
        for sample in samples
        for resolution in resolutions
    ]
    missing = [str(path) for _, _, path in rmsd_paths if not path.is_file()]
    if missing:
        raise StatisticsContractError(
            "enabled structure QC is missing required output(s): {}".format(
                ", ".join(missing)
            )
        )

    rmsd_lines = []
    for _, _, path in rmsd_paths:
        try:
            lines = path.read_text().splitlines()
        except (OSError, UnicodeError) as error:
            raise StatisticsContractError(
                "cannot read required RMSD output {!r}: {}".format(str(path), error)
            )
        matches = [line for line in lines if RMSD_RE.search(line)]
        if len(matches) != 1:
            raise StatisticsContractError(
                "required RMSD output {!r} contains {} metric lines, expected 1".format(
                    str(path), len(matches)
                )
            )
        rmsd_lines.append("{}:{}\n".format(path.relative_to(work_path), matches[0]))

    write_if_changed(output_path, "".join(sorted(rmsd_lines)).encode("utf-8"))


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_contract_arguments(command):
        command.add_argument("--contract", required=True)
        command.add_argument("--work-dir", required=True)
        command.add_argument("--output", required=True)

    manifest = subparsers.add_parser("manifest")
    add_contract_arguments(manifest)

    raw = subparsers.add_parser("raw-fastq-stat")
    add_contract_arguments(raw)
    raw.add_argument("--workers", type=int, default=1)

    rna = subparsers.add_parser("rna-alignment-stat")
    add_contract_arguments(rna)

    pairs = subparsers.add_parser("cleaned-pair-stats")
    pairs.add_argument("--contract", required=True)
    pairs.add_argument("--work-dir", required=True)
    pairs.add_argument("--pairs-output", required=True)
    pairs.add_argument("--inter-output", required=True)
    pairs.add_argument("--workers", type=int, default=1)

    rmsd = subparsers.add_parser("rmsd-stats")
    rmsd.add_argument("--contract", required=True)
    rmsd.add_argument("--work-dir", required=True)
    rmsd.add_argument("--resolutions", nargs="+", required=True)
    rmsd.add_argument("--output", required=True)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        contract = load_frozen_contract(args.contract, args.work_dir)
        if args.command == "manifest":
            write_manifest(contract, args.output)
        elif args.command == "raw-fastq-stat":
            write_raw_fastq_stats(contract, args.output, args.workers)
        elif args.command == "rna-alignment-stat":
            write_rna_alignment_stats(contract, sys.stdin, args.output)
        elif args.command == "cleaned-pair-stats":
            write_cleaned_pair_stats(
                contract,
                args.work_dir,
                args.pairs_output,
                args.inter_output,
                args.workers,
            )
        else:
            write_rmsd_stats(
                contract, args.work_dir, args.resolutions, args.output
            )
    except (InputContractError, OSError, StatisticsContractError, ValueError) as error:
        print("statistics contract error: {}".format(error), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

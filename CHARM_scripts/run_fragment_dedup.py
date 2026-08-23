#!/usr/bin/env python3
"""Run endpoint dedup atomically and report raw/final fragment counts."""

import argparse
import csv
import gzip
import hashlib
import os
import shlex
import subprocess
import sys
from pathlib import Path


def digest_and_count(path):
    digest = hashlib.sha256()
    count = 0
    with gzip.open(str(path), "rb") as handle:
        for line in handle:
            digest.update(line)
            count += 1
    return digest.hexdigest(), count


def format_command(command):
    return " ".join(shlex.quote(value) for value in command)


def write_summary(path, args, raw_digest, raw_rows, final_digest, final_rows):
    fields = [
        "sample",
        "split",
        "contract_version",
        "endpoint_eps_bp",
        "raw_rows",
        "final_deduplicated_rows",
        "raw_decompressed_sha256",
        "final_decompressed_sha256",
    ]
    row = {
        "sample": args.sample,
        "split": args.split,
        "contract_version": "r2_5prime_unclipped_v2_no_hickit_gate",
        "endpoint_eps_bp": args.eps,
        "raw_rows": raw_rows,
        "final_deduplicated_rows": final_rows,
        "raw_decompressed_sha256": raw_digest,
        "final_decompressed_sha256": final_digest,
    }
    with open(str(path), "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerow(row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--split", required=True, choices=("ct", "atac"))
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--production-script", required=True)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--eps", type=int, default=1)
    args = parser.parse_args()
    if args.eps != 1:
        parser.error("--eps must remain 1 for the frozen endpoint-dedup contract")

    input_path = Path(args.input)
    output = Path(args.output)
    stats = Path(args.stats)
    production_script = Path(args.production_script)
    log = Path(args.log)
    for required in (input_path, production_script):
        if not required.is_file():
            raise FileNotFoundError(str(required))

    output.parent.mkdir(parents=True, exist_ok=True)
    stats.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    output_tmp = Path(str(output) + ".tmp.gz")
    stats_tmp = Path(str(stats) + ".tmp")
    for temporary in (output_tmp, stats_tmp):
        if temporary.exists():
            temporary.unlink()

    raw_digest, raw_rows = digest_and_count(input_path)
    try:
        with open(str(log), "w") as log_handle:
            if raw_rows == 0:
                with gzip.open(str(output_tmp), "wb"):
                    pass
                log_handle.write("EMPTY_INPUT\t{}\t{}\n".format(args.sample, args.split))
            else:
                command = [
                    sys.executable,
                    str(production_script),
                    "-i",
                    str(input_path),
                    "-o",
                    str(output_tmp),
                    "-t",
                    "normal",
                    "-e",
                    str(args.eps),
                    "-q",
                    "0",
                ]
                log_handle.write("COMMAND\t{}\n".format(format_command(command)))
                log_handle.flush()
                subprocess.run(
                    command,
                    check=True,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                )

        final_digest, final_rows = digest_and_count(output_tmp)
        if final_rows > raw_rows:
            raise AssertionError(
                "dedup rows exceed raw rows: {} > {}".format(final_rows, raw_rows)
            )
        write_summary(
            stats_tmp,
            args,
            raw_digest,
            raw_rows,
            final_digest,
            final_rows,
        )
        os.replace(str(output_tmp), str(output))
        os.replace(str(stats_tmp), str(stats))
    finally:
        for temporary in (output_tmp, stats_tmp):
            if temporary.exists():
                temporary.unlink()

    print(
        "PASS {} {} raw_rows={} final_rows={}".format(
            args.sample, args.split, raw_rows, final_rows
        )
    )


if __name__ == "__main__":
    main()

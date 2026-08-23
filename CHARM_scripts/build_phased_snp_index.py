#!/usr/bin/env python3
"""Build a compact memory-mapped index from hickit phased-SNP TSV input."""

import argparse
from array import array
import csv
import hashlib
import os
from pathlib import Path
import re
import sys
import tempfile

from r2_fragment_selector_core import SNP_INDEX_VERSION


MANIFEST_FIELDS = [
    "index_version",
    "source_path",
    "source_size",
    "source_sha256",
    "source_data_rows",
    "source_skipped_rows",
    "chrom",
    "count",
    "positions_file",
    "alleles_file",
    "positions_sha256",
    "alleles_sha256",
]


def sha256_file(path):
    digest = hashlib.sha256()
    with open(str(path), "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_chrom_name(index, chrom):
    token = re.sub(r"[^A-Za-z0-9_.-]", "_", chrom)
    return "{:03d}.{}".format(index, token)


def write_chromosome(temp_dir, index, chrom, positions, alleles):
    prefix = safe_chrom_name(index, chrom)
    positions_name = prefix + ".positions.u32"
    alleles_name = prefix + ".alleles.u8"
    positions_path = temp_dir / positions_name
    alleles_path = temp_dir / alleles_name
    if sys.byteorder != "little":
        positions.byteswap()
    with open(str(positions_path), "wb") as handle:
        positions.tofile(handle)
    with open(str(alleles_path), "wb") as handle:
        handle.write(alleles)
    return {
        "chrom": chrom,
        "count": len(positions),
        "positions_file": positions_name,
        "alleles_file": alleles_name,
        "positions_sha256": sha256_file(positions_path),
        "alleles_sha256": sha256_file(alleles_path),
    }


def build_index(source, manifest):
    source = Path(source).resolve()
    manifest = Path(manifest)
    if not source.is_file():
        raise FileNotFoundError(str(source))
    manifest.parent.mkdir(parents=True, exist_ok=True)

    source_digest = hashlib.sha256()
    source_rows = 0
    skipped_rows = 0
    seen_chromosomes = set()
    current_chrom = None
    positions = array("I")
    alleles = bytearray()
    chrom_rows = []
    last_position = None

    with tempfile.TemporaryDirectory(
        prefix="snp_index.", dir=str(manifest.parent)
    ) as temporary_name:
        temporary = Path(temporary_name)
        with open(str(source), "rb") as handle:
            for line_number, raw_line in enumerate(handle, 1):
                source_digest.update(raw_line)
                stripped = raw_line.rstrip(b"\r\n")
                if not stripped or stripped.startswith(b"#"):
                    continue
                fields = stripped.split(b"\t")
                if len(fields) < 4:
                    raise ValueError(
                        "invalid phased-SNP row {}: expected at least 4 fields".format(
                            line_number
                        )
                    )
                source_rows += 1
                if len(fields[2]) != 1 or len(fields[3]) != 1:
                    skipped_rows += 1
                    continue
                chrom = fields[0].decode("ascii")
                position = int(fields[1]) - 1
                if position < 0 or position > 0xFFFFFFFF:
                    raise ValueError(
                        "SNP position outside uint32 range at row {}".format(line_number)
                    )
                if current_chrom != chrom:
                    if current_chrom is not None:
                        chrom_rows.append(
                            write_chromosome(
                                temporary,
                                len(chrom_rows),
                                current_chrom,
                                positions,
                                alleles,
                            )
                        )
                        seen_chromosomes.add(current_chrom)
                    if chrom in seen_chromosomes:
                        raise ValueError(
                            "chromosome {} occurs in multiple blocks".format(chrom)
                        )
                    current_chrom = chrom
                    positions = array("I")
                    alleles = bytearray()
                    last_position = None
                if last_position is not None and position < last_position:
                    raise ValueError(
                        "SNP positions are not sorted for {} at row {}".format(
                            chrom, line_number
                        )
                    )
                positions.append(position)
                alleles.extend((fields[2][0], fields[3][0]))
                last_position = position

        if current_chrom is not None:
            chrom_rows.append(
                write_chromosome(
                    temporary,
                    len(chrom_rows),
                    current_chrom,
                    positions,
                    alleles,
                )
            )
        if not chrom_rows:
            raise ValueError("no single-base phased SNPs in {}".format(source))

        common = {
            "index_version": SNP_INDEX_VERSION,
            "source_path": str(source),
            "source_size": source.stat().st_size,
            "source_sha256": source_digest.hexdigest(),
            "source_data_rows": source_rows,
            "source_skipped_rows": skipped_rows,
        }
        manifest_tmp = temporary / "manifest.tsv"
        with open(str(manifest_tmp), "w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=MANIFEST_FIELDS)
            writer.writeheader()
            for chrom_row in chrom_rows:
                row = dict(common)
                row.update(chrom_row)
                writer.writerow(row)

        for chrom_row in chrom_rows:
            for key in ("positions_file", "alleles_file"):
                os.replace(
                    str(temporary / chrom_row[key]), str(manifest.parent / chrom_row[key])
                )
        os.replace(str(manifest_tmp), str(manifest))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True, help="Output manifest.tsv")
    args = parser.parse_args()
    build_index(args.input, args.output)


if __name__ == "__main__":
    main()


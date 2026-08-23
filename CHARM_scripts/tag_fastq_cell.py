#!/usr/bin/env python3
"""Insert a safe cell code before the eight-base UMI in FASTQ read IDs.

The header transformation follows the validated logic in Part 2
``rna_read_mode/scripts/sync_tag_fastq.py``, narrowed here to the R1 FASTQ
produced directly by ``umi_tools extract``.
"""

import argparse
import gzip
import os
import re
import sys
import tempfile


SAFE_CODE_PATTERN = re.compile(br"^C[0-9a-f]{16}$")
UMI_PATTERN = re.compile(br"^[ACGTNacgtn]{8}$")


def _without_line_ending(line):
    return line.rstrip(b"\r\n")


def tagged_header(header, cell_code):
    code = cell_code.encode("ascii") if isinstance(cell_code, str) else cell_code
    if not SAFE_CODE_PATTERN.match(code):
        raise ValueError("invalid safe cell code: {!r}".format(cell_code))
    if not header.startswith(b"@"):
        raise ValueError("FASTQ header does not start with @")

    newline = b"\r\n" if header.endswith(b"\r\n") else b"\n"
    body = _without_line_ending(header)[1:]
    match = re.match(br"^(\S+)(.*)$", body)
    if match is None:
        raise ValueError("empty FASTQ read ID")
    identifier, suffix = match.groups()
    if b"_" not in identifier:
        raise ValueError(
            "UMI suffix is absent from read ID: {}".format(
                identifier.decode("ascii", errors="replace")
            )
        )
    prefix, umi = identifier.rsplit(b"_", 1)
    if not UMI_PATTERN.match(umi):
        raise ValueError(
            "invalid eight-base UMI in read ID: {}".format(
                identifier.decode("ascii", errors="replace")
            )
        )
    return b"@" + prefix + b"_" + code + b"_" + umi + suffix + newline


def tag_fastq(input_path, output_path, cell_code):
    if os.path.abspath(input_path) == os.path.abspath(output_path):
        raise ValueError("input and output FASTQ paths must differ")
    output_dir = os.path.dirname(os.path.abspath(output_path))
    os.makedirs(output_dir, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=".tag-fastq.", dir=output_dir)
    os.close(descriptor)
    records = 0
    try:
        with gzip.open(input_path, "rb") as source, open(temporary, "wb") as raw_output:
            with gzip.GzipFile(
                fileobj=raw_output, filename="", mode="wb", mtime=0
            ) as destination:
                while True:
                    header = source.readline()
                    if not header:
                        break
                    sequence = source.readline()
                    plus = source.readline()
                    quality = source.readline()
                    if not sequence or not plus or not quality:
                        raise ValueError(
                            "truncated FASTQ record {} in {}".format(
                                records + 1, input_path
                            )
                        )
                    if not plus.startswith(b"+"):
                        raise ValueError(
                            "invalid FASTQ plus line at record {} in {}".format(
                                records + 1, input_path
                            )
                        )
                    if len(_without_line_ending(sequence)) != len(
                        _without_line_ending(quality)
                    ):
                        raise ValueError(
                            "sequence/quality length mismatch at record {} in {}".format(
                                records + 1, input_path
                            )
                        )
                    destination.write(tagged_header(header, cell_code))
                    destination.write(sequence)
                    destination.write(plus)
                    destination.write(quality)
                    records += 1
            raw_output.flush()
            os.fsync(raw_output.fileno())
        os.replace(temporary, output_path)
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise
    return records


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--cell-code", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        tag_fastq(args.input, args.output, args.cell_code)
    except (OSError, UnicodeError, ValueError) as error:
        print("FASTQ cell-tagging error: {}".format(error), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Compare canonical 2D outputs while ignoring gzip bytes and record order."""

import argparse
from collections import Counter
import gzip
from pathlib import Path
import sys


ARTIFACTS = (
    "contacts.seg.gz",
    "raw.pairs.gz",
    "contacts.pairs.gz",
    "c1.pairs.gz",
    "c12.pairs.gz",
)


def _lines(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(str(path), "rt") as handle:
        return [line.rstrip("\n") for line in handle if line.rstrip("\n")]


def canonical_records(path):
    return Counter(_lines(path))


def data_qnames(path):
    return {
        line.split("\t", 1)[0]
        for line in _lines(path)
        if not line.startswith("#")
    }


def data_qname_counts(path):
    return Counter(
        line.split("\t", 1)[0]
        for line in _lines(path)
        if not line.startswith("#")
    )


def phase_records(path):
    records = Counter()
    for line in _lines(path):
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        qname = fields[0]
        for segment in fields[1:]:
            parts = segment.split("!")
            if len(parts) < 5:
                raise ValueError("malformed SEG field in {}: {}".format(path, segment))
            records[(qname, parts[0], parts[1], parts[2], parts[3], parts[4])] += 1
    return records


def compare_output_dirs(baseline, candidate, threshold=0.015):
    checks = []
    for name in ARTIFACTS:
        observed = canonical_records(candidate / name)
        expected = canonical_records(baseline / name)
        checks.append((name + ":canonical_multiset", observed == expected))
        checks.append(
            (
                name + ":qname_set",
                data_qnames(candidate / name) == data_qnames(baseline / name),
            )
        )

    checks.append(
        (
            "contacts.seg.gz:phase_records",
            phase_records(candidate / "contacts.seg.gz")
            == phase_records(baseline / "contacts.seg.gz"),
        )
    )

    for root_name, root in (("baseline", baseline), ("candidate", candidate)):
        seg_qname_counts = data_qname_counts(root / "contacts.seg.gz")
        dedup_qnames = data_qnames(root / "contacts.pairs.gz")
        raw_qnames = data_qnames(root / "raw.pairs.gz")
        c1_qnames = data_qnames(root / "c1.pairs.gz")
        c12_qnames = data_qnames(root / "c12.pairs.gz")
        checks.append(
            (
                root_name + ":one_seg_record_per_qname",
                all(count == 1 for count in seg_qname_counts.values()),
            )
        )
        checks.append((root_name + ":dedup_subset", dedup_qnames <= raw_qnames))
        checks.append((root_name + ":c1_subset", c1_qnames <= dedup_qnames))
        checks.append((root_name + ":c12_subset", c12_qnames <= c1_qnames))

    baseline_ratio = float((baseline / "yperx.txt").read_text().strip())
    candidate_ratio = float((candidate / "yperx.txt").read_text().strip())
    checks.append(("yperx:exact", baseline_ratio == candidate_ratio))
    checks.append(
        (
            "yperx:classification",
            (baseline_ratio > threshold) == (candidate_ratio > threshold),
        )
    )
    return checks


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--yperx-threshold", type=float, default=0.015)
    args = parser.parse_args(argv)

    checks = compare_output_dirs(args.baseline, args.candidate, args.yperx_threshold)
    print("check\tstatus")
    for name, passed in checks:
        print("{}\t{}".format(name, "PASS" if passed else "FAIL"))
    return 0 if all(passed for _, passed in checks) else 1


if __name__ == "__main__":
    sys.exit(main())

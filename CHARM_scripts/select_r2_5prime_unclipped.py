#!/usr/bin/env python3
"""Select physical-R2 5'-unclipped fragments with independently phased backends."""

import argparse
import csv
import gzip
import hashlib
import os
import shlex
import subprocess
import tempfile
import time
from pathlib import Path

import pysam

from r2_fragment_selector_core import (
    CONTRACT_VERSION,
    PhasedSnpIndex,
    annotate_record_geometry,
    annotation_priority,
    bed_line,
    iter_groups,
    phase_record_bounds,
    select_fragment,
    set_annotation_phase,
)


ANNOTATION_FIELDS = [
    "qname",
    "origin",
    "flag",
    "chrom",
    "start",
    "end",
    "strand",
    "mapq",
    "q_start",
    "q_end",
    "aligned_query_length",
    "phase",
    "phase0_count",
    "phase1_count",
    "other_allele_count",
]
INTEGER_ANNOTATION_FIELDS = {
    "flag",
    "start",
    "end",
    "mapq",
    "q_start",
    "q_end",
    "aligned_query_length",
    "phase0_count",
    "phase1_count",
    "other_allele_count",
}
PHASE_FIELDS = [
    "sample",
    "split",
    "qname",
    "chrom",
    "start",
    "end",
    "strand",
    "mapq",
    "q_start",
    "selected_flag",
    "selected_alignment_type",
    "r2_candidate_count",
    "direct_r2_phase",
    "r1_fallback_phase",
    "final_phase",
    "phase_source",
    "r2_phase0_snp_count",
    "r2_phase1_snp_count",
    "r2_other_allele_count",
    "informative_r1_alignment_count",
    "r1_internal_conflict",
    "r2_r1_conflict",
]
SUMMARY_FIELDS = [
    "sample",
    "split",
    "contract_version",
    "min_mapq",
    "min_baseq",
    "bam_qnames",
    "qnames_with_any_mapq20_alignment",
    "annotated_alignments",
    "qnames_with_mapq20_nonsecondary_r2",
    "selected_first_r2",
    "selected_first_r2_clipped",
    "r2_5prime_unclipped_raw_rows",
    "direct_r2_phased_rows",
    "r1_fallback_phased_rows",
    "unphased_rows",
    "r1_internal_conflict_rows",
    "r2_r1_conflict_rows",
    "raw_decompressed_sha256",
    "phase_decompressed_sha256",
]


def format_command(command):
    return " ".join(shlex.quote(str(value)) for value in command)


def parse_annotation(fields, line_number):
    if len(fields) != len(ANNOTATION_FIELDS):
        raise ValueError(
            "sam2phase line {} has {} fields; expected {}".format(
                line_number, len(fields), len(ANNOTATION_FIELDS)
            )
        )
    annotation = dict(zip(ANNOTATION_FIELDS, fields))
    for field in INTEGER_ANNOTATION_FIELDS:
        annotation[field] = int(annotation[field])
    if annotation["phase"] not in (".", "0", "1"):
        raise ValueError(
            "invalid phase {!r} at sam2phase line {}".format(
                annotation["phase"], line_number
            )
        )
    return annotation


def iter_reference_annotations(path):
    header_seen = False
    with open(str(path)) as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#columns:"):
                observed = line[len("#columns:") :].strip().split("\t")
                if observed != ANNOTATION_FIELDS:
                    raise ValueError(
                        "unexpected sam2phase columns: {}".format(" ".join(observed))
                    )
                header_seen = True
                continue
            if line.startswith("#"):
                continue
            if not header_seen:
                raise ValueError("sam2phase data precede #columns header")
            yield parse_annotation(line.split("\t"), line_number)
    if not header_seen:
        raise ValueError("sam2phase did not emit a #columns header")


def iter_annotation_groups(annotations):
    current_name = None
    group = []
    for annotation in annotations:
        qname = annotation["qname"]
        if current_name is None:
            current_name = qname
        if qname != current_name:
            yield current_name, group
            current_name = qname
            group = []
        group.append(annotation)
    if current_name is not None:
        yield current_name, group


def run_name_sort(bam, namesort_bam, threads, temporary, log_handle):
    command = [
        "samtools",
        "sort",
        "-n",
        "-@",
        str(threads),
        "-T",
        str(temporary / "sorttmp"),
        "-o",
        str(namesort_bam),
        str(bam),
    ]
    log_handle.write("COMMAND\t{}\n".format(format_command(command)))
    log_handle.flush()
    subprocess.run(command, check=True, stdout=log_handle, stderr=log_handle)


def run_hickit_reference(namesort_bam, annotation_path, hickit, snp, args, log_handle):
    env = os.environ.copy()
    env["PATH"] = str(hickit.parent) + os.pathsep + env.get("PATH", "")
    view_command = ["samtools", "view", str(namesort_bam)]
    hickit_command = [
        str(hickit),
        "sam2phase",
        "-q",
        str(args.min_mapq),
        "-Q",
        str(args.min_baseq),
        "-v",
        str(snp),
        "-",
    ]
    log_handle.write(
        "COMMAND\t{} | {}\n".format(
            format_command(view_command), format_command(hickit_command)
        )
    )
    log_handle.flush()
    with open(str(annotation_path), "wb") as output_handle:
        view = subprocess.Popen(
            view_command, stdout=subprocess.PIPE, stderr=log_handle, env=env
        )
        hickit_process = subprocess.Popen(
            hickit_command,
            stdin=view.stdout,
            stdout=output_handle,
            stderr=log_handle,
            env=env,
        )
        view.stdout.close()
        hickit_code = hickit_process.wait()
        view_code = view.wait()
    if (view_code, hickit_code) != (0, 0):
        raise RuntimeError(
            "samtools view / hickit sam2phase failed: {} {}".format(
                view_code, hickit_code
            )
        )


def count_bam_qnames(namesort_bam):
    count = 0
    with pysam.AlignmentFile(str(namesort_bam), "rb") as handle:
        for _, _ in iter_groups(handle.fetch(until_eof=True)):
            count += 1
    return count


def phase_python_batch(batch, snp_index, min_baseq):
    targets_by_chromosome = {}
    for _, entries in batch:
        r2_candidates = [
            entry
            for entry in entries
            if entry[1]["origin"] == "R2" and not entry[1]["flag"] & 0x100
        ]
        if not r2_candidates:
            continue
        selected = min(r2_candidates, key=lambda entry: annotation_priority(entry[1]))
        if selected[1]["q_start"] != 0:
            continue

        phase_targets = [selected]
        phase_targets.extend(
            entry
            for entry in entries
            if entry[1]["origin"] == "R1" and not entry[1]["flag"] & 0x100
        )
        for target in phase_targets:
            annotation = target[1]
            targets_by_chromosome.setdefault(annotation["chrom"], []).append(target)

    for chrom, targets in targets_by_chromosome.items():
        bounds = snp_index.batch_bounds(
            chrom,
            [annotation["start"] for _, annotation in targets],
            [annotation["end"] for _, annotation in targets],
        )
        if bounds is None:
            continue
        positions, alleles, lefts, rights = bounds
        for (record, annotation), left, right in zip(targets, lefts, rights):
            set_annotation_phase(
                annotation,
                phase_record_bounds(
                    record,
                    chrom,
                    positions,
                    alleles,
                    left,
                    right,
                    min_baseq,
                ),
            )

    return [
        (qname, [annotation for _, annotation in entries])
        for qname, entries in batch
    ]


def python_annotation_groups(
    namesort_bam,
    snp_index,
    min_mapq,
    min_baseq,
    batch_qnames,
    counters=None,
):
    with pysam.AlignmentFile(str(namesort_bam), "rb") as handle:
        batch = []
        for qname, records in iter_groups(handle.fetch(until_eof=True)):
            if counters is not None:
                counters["bam_qnames"] += 1
            entries = [
                (record, annotate_record_geometry(record, handle.header))
                for record in records
                if not record.is_unmapped and record.mapping_quality >= min_mapq
            ]
            batch.append((qname, entries))
            if len(batch) >= batch_qnames:
                yield from phase_python_batch(batch, snp_index, min_baseq)
                batch = []
        if batch:
            yield from phase_python_batch(batch, snp_index, min_baseq)


def phase_row(args, qname, result):
    selected = result["selected"]
    if selected["flag"] & 0x800:
        alignment_type = "supplementary"
    else:
        alignment_type = "primary"
    return {
        "sample": args.sample,
        "split": args.split,
        "qname": qname,
        "chrom": selected["chrom"],
        "start": selected["start"],
        "end": selected["end"],
        "strand": selected["strand"],
        "mapq": selected["mapq"],
        "q_start": selected["q_start"],
        "selected_flag": selected["flag"],
        "selected_alignment_type": alignment_type,
        "r2_candidate_count": result["r2_candidate_count"],
        "direct_r2_phase": result["direct_phase"],
        "r1_fallback_phase": result["r1_phase"],
        "final_phase": result["final_phase"],
        "phase_source": result["phase_source"],
        "r2_phase0_snp_count": selected["phase0_count"],
        "r2_phase1_snp_count": selected["phase1_count"],
        "r2_other_allele_count": selected["other_allele_count"],
        "informative_r1_alignment_count": result["informative_r1_count"],
        "r1_internal_conflict": int(result["r1_internal_conflict"]),
        "r2_r1_conflict": int(result["r2_r1_conflict"]),
    }


def process_groups(groups, args, output_path, phase_path, counters):
    raw_digest = hashlib.sha256()
    phase_digest = hashlib.sha256()
    with gzip.open(str(output_path), "wt", compresslevel=6) as output_handle:
        with gzip.open(str(phase_path), "wt", compresslevel=6, newline="") as phase_handle:
            phase_header = "\t".join(PHASE_FIELDS) + "\n"
            phase_handle.write(phase_header)
            phase_digest.update(phase_header.encode("utf-8"))
            for qname, annotations in groups:
                if not annotations:
                    continue
                counters["qnames_with_any_mapq20_alignment"] += 1
                counters["annotated_alignments"] += len(annotations)
                if any(
                    item["origin"] == "R2" and not item["flag"] & 0x100
                    for item in annotations
                ):
                    counters["qnames_with_mapq20_nonsecondary_r2"] += 1
                result = select_fragment(annotations)
                if result is None:
                    continue
                counters["selected_first_r2"] += 1
                selected = result["selected"]
                if selected["q_start"] != 0:
                    counters["selected_first_r2_clipped"] += 1
                    continue

                raw_line = bed_line(selected, qname)
                output_handle.write(raw_line)
                raw_digest.update(raw_line.encode("utf-8"))
                counters["r2_5prime_unclipped_raw_rows"] += 1

                row = phase_row(args, qname, result)
                phase_line = "\t".join(str(row[field]) for field in PHASE_FIELDS) + "\n"
                phase_handle.write(phase_line)
                phase_digest.update(phase_line.encode("utf-8"))
                if result["phase_source"] == "R2":
                    counters["direct_r2_phased_rows"] += 1
                elif result["phase_source"] == "R1":
                    counters["r1_fallback_phased_rows"] += 1
                else:
                    counters["unphased_rows"] += 1
                counters["r1_internal_conflict_rows"] += int(
                    result["r1_internal_conflict"]
                )
                counters["r2_r1_conflict_rows"] += int(result["r2_r1_conflict"])
    return raw_digest.hexdigest(), phase_digest.hexdigest()


def write_summary(path, args, counters, raw_sha256, phase_sha256):
    row = {
        "sample": args.sample,
        "split": args.split,
        "contract_version": CONTRACT_VERSION,
        "min_mapq": args.min_mapq,
        "min_baseq": args.min_baseq,
        "raw_decompressed_sha256": raw_sha256,
        "phase_decompressed_sha256": phase_sha256,
    }
    row.update(counters)
    with open(str(path), "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerow(row)


def build(args, log_handle):
    total_started = time.perf_counter()
    bam = Path(args.bam)
    output = Path(args.output)
    phase_output = Path(args.phase_output)
    stats = Path(args.stats)
    temp_root = Path(args.temp_root)
    if not bam.is_file():
        raise FileNotFoundError(str(bam))
    if bam.stat().st_size == 0:
        raise ValueError("input BAM is zero bytes: {}".format(bam))
    if args.backend == "python" and not args.snp_index:
        raise ValueError("--snp-index is required for --backend python")
    if args.backend == "hickit_reference" and (not args.hickit or not args.snp):
        raise ValueError(
            "--hickit and --snp are required for --backend hickit_reference"
        )

    for path in (output, phase_output, stats):
        path.parent.mkdir(parents=True, exist_ok=True)
    temp_root.mkdir(parents=True, exist_ok=True)
    temporary_outputs = [Path(str(path) + ".tmp") for path in (output, phase_output, stats)]
    for path in temporary_outputs:
        if path.exists():
            path.unlink()

    counters = {
        "bam_qnames": 0,
        "qnames_with_any_mapq20_alignment": 0,
        "annotated_alignments": 0,
        "qnames_with_mapq20_nonsecondary_r2": 0,
        "selected_first_r2": 0,
        "selected_first_r2_clipped": 0,
        "r2_5prime_unclipped_raw_rows": 0,
        "direct_r2_phased_rows": 0,
        "r1_fallback_phased_rows": 0,
        "unphased_rows": 0,
        "r1_internal_conflict_rows": 0,
        "r2_r1_conflict_rows": 0,
    }
    try:
        with tempfile.TemporaryDirectory(
            prefix="{}.{}.".format(args.sample, args.split), dir=str(temp_root)
        ) as temporary_name:
            temporary = Path(temporary_name)
            namesort_bam = temporary / "input.namesort.bam"
            stage_started = time.perf_counter()
            run_name_sort(bam, namesort_bam, args.threads, temporary, log_handle)
            name_sort_seconds = time.perf_counter() - stage_started
            qname_count_seconds = 0.0
            index_load_seconds = 0.0

            if args.backend == "python":
                stage_started = time.perf_counter()
                snp_index = PhasedSnpIndex(args.snp_index)
                index_load_seconds = time.perf_counter() - stage_started
                groups = python_annotation_groups(
                    namesort_bam,
                    snp_index,
                    args.min_mapq,
                    args.min_baseq,
                    args.python_batch_qnames,
                    counters,
                )
            else:
                stage_started = time.perf_counter()
                counters["bam_qnames"] = count_bam_qnames(namesort_bam)
                qname_count_seconds = time.perf_counter() - stage_started
                hickit = Path(args.hickit)
                snp = Path(args.snp)
                for required in (hickit, snp):
                    if not required.is_file():
                        raise FileNotFoundError(str(required))
                annotation_path = temporary / "hickit.sam2phase.tsv"
                run_hickit_reference(
                    namesort_bam,
                    annotation_path,
                    hickit,
                    snp,
                    args,
                    log_handle,
                )
                groups = iter_annotation_groups(iter_reference_annotations(annotation_path))

            stage_started = time.perf_counter()
            raw_sha256, phase_sha256 = process_groups(
                groups, args, temporary_outputs[0], temporary_outputs[1], counters
            )
            select_write_seconds = time.perf_counter() - stage_started
        write_summary(
            temporary_outputs[2], args, counters, raw_sha256, phase_sha256
        )
        for temporary_output, final_output in zip(
            temporary_outputs, (output, phase_output, stats)
        ):
            os.replace(str(temporary_output), str(final_output))
        total_seconds = time.perf_counter() - total_started
    finally:
        for temporary_output in temporary_outputs:
            if temporary_output.exists():
                temporary_output.unlink()

    for stage, seconds in (
        ("name_sort_seconds", name_sort_seconds),
        ("qname_count_seconds", qname_count_seconds),
        ("index_load_seconds", index_load_seconds),
        ("select_write_seconds", select_write_seconds),
        ("total_seconds", total_seconds),
    ):
        log_handle.write("TIMING\t{}\t{:.6f}\n".format(stage, seconds))
    log_handle.write(
        "PASS\t{}\t{}\tbackend={}\trows={}\traw_sha256={}\tphase_sha256={}\n".format(
            args.sample,
            args.split,
            args.backend,
            counters["r2_5prime_unclipped_raw_rows"],
            raw_sha256,
            phase_sha256,
        )
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--split", required=True, choices=("ct", "atac"))
    parser.add_argument("--bam", required=True)
    parser.add_argument(
        "--backend", choices=("python", "hickit_reference"), default="python"
    )
    parser.add_argument("--snp-index")
    parser.add_argument("--hickit")
    parser.add_argument("--snp")
    parser.add_argument("--output", required=True)
    parser.add_argument("--phase-output", required=True)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--temp-root", required=True)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--python-batch-qnames", type=int, default=32768)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    if args.min_mapq < 0 or args.min_baseq < 0:
        parser.error("quality thresholds must be non-negative")
    if args.python_batch_qnames < 1:
        parser.error("--python-batch-qnames must be at least 1")

    log_path = Path(args.log)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(str(log_path), "w") as log_handle:
        log_handle.write("CONTRACT\t{}\n".format(CONTRACT_VERSION))
        log_handle.write("BACKEND\t{}\n".format(args.backend))
        build(args, log_handle)


if __name__ == "__main__":
    main()

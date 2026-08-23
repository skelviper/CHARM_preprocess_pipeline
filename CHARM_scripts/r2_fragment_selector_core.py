#!/usr/bin/env python3
"""Shared contract logic for physical-R2 fragment selection and phasing."""

import csv
from pathlib import Path

import numpy as np


CONTRACT_VERSION = "r2_5prime_unclipped_v2_no_hickit_gate"
SNP_INDEX_VERSION = "phased_snp_u32_u8_v1"
SUPPORTED_CIGAR_OPERATIONS = {0, 1, 2, 4, 5}


def iter_groups(records):
    current_name = None
    group = []
    for record in records:
        if current_name is None:
            current_name = record.query_name
        if record.query_name != current_name:
            yield current_name, group
            current_name = record.query_name
            group = []
        group.append(record)
    if current_name is not None:
        yield current_name, group


def cigar_metrics(cigartuples, qname="unknown"):
    cigar = cigartuples or []
    unsupported = sorted({operation for operation, _ in cigar} - SUPPORTED_CIGAR_OPERATIONS)
    if unsupported:
        raise ValueError(
            "unsupported CIGAR operation(s) {} for {}".format(unsupported, qname)
        )

    left_clip = 0
    for operation, length in cigar:
        if operation not in (4, 5):
            break
        left_clip += length

    right_clip = 0
    for operation, length in reversed(cigar):
        if operation not in (4, 5):
            break
        right_clip += length

    reference_length = sum(
        length for operation, length in cigar if operation in (0, 2)
    )
    aligned_query_length = sum(
        length for operation, length in cigar if operation in (0, 1)
    )
    full_query_length = sum(
        length for operation, length in cigar if operation in (0, 1, 4, 5)
    )
    return {
        "left_clip": left_clip,
        "right_clip": right_clip,
        "reference_length": reference_length,
        "aligned_query_length": aligned_query_length,
        "full_query_length": full_query_length,
    }


class PhasedSnpIndex:
    """Memory-mapped phased SNP index produced by build_phased_snp_index.py."""

    def __init__(self, manifest_path):
        self.manifest_path = Path(manifest_path)
        self.chromosomes = {}
        with open(str(self.manifest_path), newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows:
            raise ValueError("empty SNP index manifest: {}".format(self.manifest_path))

        for row in rows:
            if row["index_version"] != SNP_INDEX_VERSION:
                raise ValueError(
                    "unsupported SNP index version {} in {}".format(
                        row["index_version"], self.manifest_path
                    )
                )
            chrom = row["chrom"]
            count = int(row["count"])
            positions_path = self.manifest_path.parent / row["positions_file"]
            alleles_path = self.manifest_path.parent / row["alleles_file"]
            if positions_path.stat().st_size != count * 4:
                raise ValueError("invalid SNP positions size: {}".format(positions_path))
            if alleles_path.stat().st_size != count * 2:
                raise ValueError("invalid SNP alleles size: {}".format(alleles_path))
            positions = np.memmap(
                str(positions_path), dtype="<u4", mode="r", shape=(count,)
            )
            alleles = np.memmap(
                str(alleles_path), dtype="u1", mode="r", shape=(count, 2)
            )
            self.chromosomes[chrom] = (positions, alleles)

    def overlap(self, chrom, start, end):
        indexed = self.chromosomes.get(chrom)
        if indexed is None or start >= end:
            return None, None
        positions, alleles = indexed
        left = int(np.searchsorted(positions, start, side="left"))
        right = int(np.searchsorted(positions, end, side="left"))
        return positions[left:right], alleles[left:right]

    def batch_bounds(self, chrom, starts, ends):
        """Return vectorized SNP-index bounds for one chromosome."""
        indexed = self.chromosomes.get(chrom)
        if indexed is None:
            return None
        positions, alleles = indexed
        starts_array = np.asarray(starts, dtype=np.int64)
        ends_array = np.asarray(ends, dtype=np.int64)
        lefts = np.searchsorted(positions, starts_array, side="left")
        rights = np.searchsorted(positions, ends_array, side="left")
        return positions, alleles, lefts, rights


def phase_record_bounds(
    record, chrom, positions, alleles, left, right, min_baseq
):
    """Phase one record using SNP bounds already computed for its full span."""
    left = int(left)
    right = int(right)
    if left >= right:
        return ".", 0, 0, 0

    sequence = record.query_sequence or ""
    qualities = record.query_qualities
    use_qualities = qualities is not None and len(qualities) == len(sequence)
    reference_position = record.reference_start
    query_position = 0
    variant_index = left
    phase0_count = 0
    phase1_count = 0
    other_allele_count = 0

    for operation, length in record.cigartuples or []:
        if operation == 0:
            reference_end = reference_position + length
            while (
                variant_index < right
                and int(positions[variant_index]) < reference_position
            ):
                variant_index += 1
            while variant_index < right:
                position = int(positions[variant_index])
                if position >= reference_end:
                    break
                query_index = query_position + position - reference_position
                if query_index < 0 or query_index >= len(sequence):
                    raise ValueError(
                        "CIGAR parsing error for {} at {}:{}".format(
                            record.query_name, chrom, position + 1
                        )
                    )
                base_quality = qualities[query_index] if use_qualities else min_baseq
                if base_quality >= min_baseq:
                    base = ord(sequence[query_index])
                    if base == int(alleles[variant_index, 0]):
                        phase0_count += 1
                    elif base == int(alleles[variant_index, 1]):
                        phase1_count += 1
                    else:
                        other_allele_count += 1
                variant_index += 1
            reference_position = reference_end
            query_position += length
        elif operation == 1:
            query_position += length
        elif operation == 2:
            reference_position += length
        elif operation == 4:
            query_position += length
        elif operation == 5:
            pass
        else:
            raise ValueError(
                "unsupported CIGAR operation {} for {}".format(
                    operation, record.query_name
                )
            )

    if phase0_count and not phase1_count:
        phase = "0"
    elif phase1_count and not phase0_count:
        phase = "1"
    else:
        phase = "."
    return phase, phase0_count, phase1_count, other_allele_count


def phase_record(record, chrom, snp_index, min_baseq):
    indexed = snp_index.chromosomes.get(chrom)
    if indexed is None:
        return ".", 0, 0, 0
    positions, alleles = indexed
    left = np.searchsorted(positions, record.reference_start, side="left")
    right = np.searchsorted(positions, record.reference_end, side="left")
    return phase_record_bounds(
        record, chrom, positions, alleles, left, right, min_baseq
    )


def annotate_record_geometry(record, header):
    metrics = cigar_metrics(record.cigartuples, record.query_name)
    chrom = header.get_reference_name(record.reference_id)
    origin = "R1" if record.is_read1 else "R2" if record.is_read2 else "."
    q_start = metrics["right_clip"] if record.is_reverse else metrics["left_clip"]
    q_end = (
        metrics["full_query_length"] - metrics["left_clip"]
        if record.is_reverse
        else metrics["full_query_length"] - metrics["right_clip"]
    )

    if origin == "R2":
        strand = "+" if record.is_reverse else "-"
    else:
        strand = "-" if record.is_reverse else "+"

    return {
        "qname": record.query_name,
        "origin": origin,
        "flag": record.flag,
        "chrom": chrom,
        "start": record.reference_start,
        "end": record.reference_start + metrics["reference_length"],
        "strand": strand,
        "mapq": record.mapping_quality,
        "q_start": q_start,
        "q_end": q_end,
        "aligned_query_length": metrics["aligned_query_length"],
        "phase": ".",
        "phase0_count": 0,
        "phase1_count": 0,
        "other_allele_count": 0,
    }


def set_annotation_phase(annotation, phase_result):
    (
        annotation["phase"],
        annotation["phase0_count"],
        annotation["phase1_count"],
        annotation["other_allele_count"],
    ) = phase_result


def annotate_record(record, header, snp_index, min_baseq):
    annotation = annotate_record_geometry(record, header)
    set_annotation_phase(
        annotation,
        phase_record(record, annotation["chrom"], snp_index, min_baseq),
    )
    return annotation


def annotation_priority(annotation):
    return (
        annotation["q_start"],
        1 if annotation["flag"] & 0x800 else 0,
        -annotation["mapq"],
        -annotation["aligned_query_length"],
        annotation["chrom"],
        annotation["start"],
        annotation["end"],
        annotation["strand"],
        annotation["flag"],
    )


def select_fragment(annotations):
    r2_candidates = [
        annotation
        for annotation in annotations
        if annotation["origin"] == "R2" and not annotation["flag"] & 0x100
    ]
    if not r2_candidates:
        return None

    selected = min(r2_candidates, key=annotation_priority)
    known_r1 = sorted(
        (
            annotation
            for annotation in annotations
            if annotation["origin"] == "R1"
            and not annotation["flag"] & 0x100
            and annotation["phase"] in ("0", "1")
        ),
        key=annotation_priority,
    )
    r1_phase = known_r1[0]["phase"] if known_r1 else "."
    r1_internal_conflict = len({item["phase"] for item in known_r1}) > 1
    direct_phase = selected["phase"]
    r2_r1_conflict = (
        direct_phase in ("0", "1")
        and r1_phase in ("0", "1")
        and direct_phase != r1_phase
    )
    if direct_phase in ("0", "1"):
        final_phase = direct_phase
        phase_source = "R2"
    elif r1_phase in ("0", "1"):
        final_phase = r1_phase
        phase_source = "R1"
    else:
        final_phase = "."
        phase_source = "unphased"

    return {
        "selected": selected,
        "r2_candidate_count": len(r2_candidates),
        "r1_phase": r1_phase,
        "informative_r1_count": len(known_r1),
        "r1_internal_conflict": r1_internal_conflict,
        "direct_phase": direct_phase,
        "final_phase": final_phase,
        "phase_source": phase_source,
        "r2_r1_conflict": r2_r1_conflict,
    }


def bed_line(selected, qname):
    return "{}\t{}\t{}\t{}\t{}\t{}\n".format(
        selected["chrom"],
        selected["start"],
        selected["end"],
        qname,
        selected["mapq"],
        selected["strand"],
    )

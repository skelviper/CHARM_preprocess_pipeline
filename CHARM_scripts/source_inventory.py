#!/usr/bin/env python3
"""Deterministic inventory of the executable CHARM pipeline source surface."""

import hashlib
import os
import tempfile
from collections import OrderedDict


class SourceInventoryError(RuntimeError):
    pass


FIXED_SOURCE_GROUPS = OrderedDict(
    [
        (
            "pipeline",
            (
                "runCHARM.sh",
                "CHARM.smk",
                "config.yaml",
            ),
        ),
        (
            "rules",
            (
                "rules/CHARM_split.rules",
                "rules/CHARM_cuttag.rules",
                "rules/scHiC_2dprocess.rules",
                "rules/scHiC_3dprocess.rules",
                "rules/CHARM_RNA.rules",
            ),
        ),
        (
            "helpers",
            (
                "CHARM_scripts/CHARM_dedup.R",
                "CHARM_scripts/CTHiRES.extract_dedup_reads.R",
                "CHARM_scripts/atac.smk",
                "CHARM_scripts/aggregate_inputs.py",
                "CHARM_scripts/audit_complete_run.py",
                "CHARM_scripts/build_phased_snp_index.py",
                "CHARM_scripts/classify_2d_bam.py",
                "CHARM_scripts/classify_rna_gene_compatibility.py",
                "CHARM_scripts/config.yaml",
                "CHARM_scripts/dedup_bed.py",
                "CHARM_scripts/detect_adapter.py",
                "CHARM_scripts/fanout_cutadapt.py",
                "CHARM_scripts/fraction.R",
                "CHARM_scripts/generateColor2.py",
                "CHARM_scripts/generateStat.sh",
                "CHARM_scripts/generate_stat_contract.py",
                "CHARM_scripts/hickit",
                "CHARM_scripts/hickit.js",
                "CHARM_scripts/hires3dAligner.py",
                "CHARM_scripts/input_contract.py",
                "CHARM_scripts/k8",
                "CHARM_scripts/metadata_qc.R",
                "CHARM_scripts/output_utils.sh",
                "CHARM_scripts/r2_fragment_selector_core.py",
                "CHARM_scripts/rename_count_matrix.py",
                "CHARM_scripts/resolve_run_config.py",
                "CHARM_scripts/run_fragment_dedup.py",
                "CHARM_scripts/run_cutadapt46.py",
                "CHARM_scripts/seg2bed.py",
                "CHARM_scripts/select_r2_5prime_unclipped.py",
                "CHARM_scripts/source_inventory.py",
                "CHARM_scripts/stat.ipynb",
                "CHARM_scripts/summarize_rna_output_modes.py",
                "CHARM_scripts/snakemake_compat.py",
                "CHARM_scripts/tag_fastq_cell.py",
                "CHARM_scripts/tss.py",
            ),
        ),
        (
            "environment",
            (
                "envs/charm.yml",
            ),
        ),
    ]
)

CHARMTOOLS_EXECUTION_MEMBERS = (
    "__init__.py",
    "__main__.py",
    "cli.py",
    "charm_preprocess/clean3.py",
    "charm_preprocess/clean_isolated.py",
    "charm_preprocess/clean_leg.py",
    "charm_preprocess/clean_splicing.py",
    "utils/CHARMio.py",
    "ref/__init__.py",
    "ref/chrom_alias.csv",
)


def _charmtools_manifest_members(manifest):
    if not os.path.isfile(manifest):
        return []
    members = []
    try:
        with open(manifest, encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, 1):
                line = raw_line.rstrip("\n")
                if not line:
                    continue
                fields = line.split(None, 1)
                if len(fields) != 2 or len(fields[0]) != 64:
                    raise SourceInventoryError(
                        "invalid CHARMtools manifest row {}:{}".format(
                            manifest, line_number
                        )
                    )
                relative = fields[1].lstrip("*")
                normalized = os.path.normpath(relative)
                if (
                    os.path.isabs(relative)
                    or normalized in ("", ".", "..")
                    or normalized.startswith(".." + os.sep)
                ):
                    raise SourceInventoryError(
                        "unsafe CHARMtools manifest path {}:{}: {}".format(
                            manifest, line_number, relative
                        )
                    )
                members.append(normalized)
    except (OSError, UnicodeError, SourceInventoryError):
        # A missing sentinel becomes a target-scoped MissingInput only when
        # provenance is requested; unrelated module targets can still build.
        return [".invalid-UPSTREAM_MANIFEST.sha256"]
    return sorted(set(members))


def _add_charmtools_runtime(by_path, charmtools_root, category):
    manifest = os.path.join(charmtools_root, "UPSTREAM_MANIFEST.sha256")
    by_path.setdefault(os.path.abspath(manifest), set()).add("charmtools_identity")
    provenance = os.path.join(charmtools_root, "UPSTREAM_PROVENANCE.md")
    by_path.setdefault(os.path.abspath(provenance), set()).add("charmtools_identity")
    manifest_members = set(_charmtools_manifest_members(manifest))
    if not set(CHARMTOOLS_EXECUTION_MEMBERS).issubset(manifest_members):
        invalid = os.path.join(charmtools_root, ".invalid-UPSTREAM_MANIFEST.sha256")
        by_path.setdefault(os.path.abspath(invalid), set()).add("charmtools_identity")
    for relative in CHARMTOOLS_EXECUTION_MEMBERS:
        absolute = os.path.normpath(os.path.join(charmtools_root, relative))
        runtime_category = (
            "charmtools_runtime_data"
            if relative == "ref/chrom_alias.csv"
            else category
        )
        by_path.setdefault(absolute, set()).add(runtime_category)


def build_source_inventory(pipeline_dir, config_sources=(), charmtools_dir=None):
    pipeline_dir = os.path.abspath(pipeline_dir)
    by_path = {}
    for category, relative_paths in FIXED_SOURCE_GROUPS.items():
        for relative in relative_paths:
            absolute = os.path.normpath(os.path.join(pipeline_dir, relative))
            by_path.setdefault(absolute, set()).add(category)

    if charmtools_dir is None:
        charmtools_dir = os.path.join(os.path.dirname(pipeline_dir), "CHARMtools")
    _add_charmtools_runtime(
        by_path, os.path.abspath(charmtools_dir), "charmtools_runtime"
    )

    for config_source in config_sources:
        absolute = os.path.abspath(config_source)
        by_path.setdefault(absolute, set()).add("effective_config_source")

    return [
        (",".join(sorted(categories)), path)
        for path, categories in sorted(by_path.items(), key=lambda item: item[0])
    ]


def _sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_text_if_changed(output, content):
    output = os.path.abspath(output)
    parent = os.path.dirname(output)
    if parent:
        os.makedirs(parent, exist_ok=True)
    if os.path.exists(output):
        with open(output, encoding="utf-8") as handle:
            if handle.read() == content:
                return False
    descriptor, temporary = tempfile.mkstemp(
        prefix=os.path.basename(output) + ".", dir=parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
        os.replace(temporary, output)
        return True
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def write_source_hash_inventory(entries, pipeline_dir, output):
    pipeline_dir = os.path.abspath(pipeline_dir)
    rendered = []
    for category, path in entries:
        if not os.path.exists(path):
            raise SourceInventoryError("expected pipeline source is missing: {}".format(path))
        if not os.path.isfile(path):
            raise SourceInventoryError(
                "expected pipeline source is not a regular file: {}".format(path)
            )
        if not os.access(path, os.R_OK):
            raise SourceInventoryError("expected pipeline source is unreadable: {}".format(path))
        relative = os.path.relpath(path, pipeline_dir)
        rendered.append(
            (
                relative,
                "{}\t{}\t{}".format(category, _sha256_file(path), relative),
            )
        )
    rows = ["category\tsha256\tpath"] + [
        row for _, row in sorted(rendered, key=lambda item: item[0])
    ]
    return _write_text_if_changed(output, "\n".join(rows) + "\n")

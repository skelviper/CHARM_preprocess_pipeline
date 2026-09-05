####################################
#       CHARM_pipline              #
#@author Z Liu                     #
#@Ver 0.3.0                        #
#@date 2023/8/7                    #
####################################

#############CONFIG#################

import hashlib
import json
import os
import sys

PIPELINE_DIR = workflow.basedir
SCRIPTS_DIR = os.path.join(PIPELINE_DIR, "CHARM_scripts")

if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from input_contract import InputContractError, load_and_validate_contract

CONFIG_PATH = os.path.join(PIPELINE_DIR, "config.yaml")
configfile: CONFIG_PATH

EXPERIMENT_TYPE = str(config.get("experiment_type", "")).strip().lower()
if EXPERIMENT_TYPE not in {"charm", "hires"}:
    raise ValueError(
        "experiment_type must be exactly 'charm' or 'hires', observed {!r}".format(
            config.get("experiment_type")
        )
    )
config["experiment_type"] = EXPERIMENT_TYPE
IS_CHARM = EXPERIMENT_TYPE == "charm"

CONFIG_SOURCE_FILES = list(
    dict.fromkeys(
        os.path.abspath(path)
        for path in [CONFIG_PATH] + list(workflow.overwrite_configfiles or [])
    )
)

CHARMTOOLS_DIR = config.get("softwares", {}).get(
    "CHARMtools", os.path.join(PIPELINE_DIR, "..", "CHARMtools")
)
if not os.path.isabs(CHARMTOOLS_DIR):
    CHARMTOOLS_DIR = os.path.normpath(os.path.join(PIPELINE_DIR, CHARMTOOLS_DIR))
CHARMTOOLS_PYTHONPATH = os.path.dirname(CHARMTOOLS_DIR)

if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)
from source_inventory import build_source_inventory, write_source_hash_inventory


LAUNCHER_WORK_DIR = os.environ.get("CHARM_EFFECTIVE_WORK_DIR")
if LAUNCHER_WORK_DIR:
    if not os.path.isabs(LAUNCHER_WORK_DIR):
        raise ValueError("CHARM_EFFECTIVE_WORK_DIR must be absolute")
    WORK_DIR = os.path.normpath(LAUNCHER_WORK_DIR)
else:
    WORK_DIR = config.get("work_dir", "..")
    if isinstance(WORK_DIR, str):
        WORK_DIR = WORK_DIR.replace("\r", "").strip()
    if not os.path.isabs(WORK_DIR):
        WORK_DIR = os.path.normpath(os.path.join(PIPELINE_DIR, WORK_DIR))
workdir: WORK_DIR
config["work_dir"] = WORK_DIR

# Input discovery is centralized in the launcher-generated contract. Included
# rules consume only these frozen paths and never guess mates independently.
INPUT_CONTRACT_PATH = os.path.join(
    WORK_DIR, "qc", "input_contract", "current.json"
)
try:
    INPUT_CONTRACT = load_and_validate_contract(WORK_DIR)
except InputContractError as error:
    raise ValueError(
        "invalid or missing Rawdata input contract: {}. Run runCHARM.sh first.".format(
            error
        )
    )

SAMPLES = [sample["sample_name"] for sample in INPUT_CONTRACT["samples"]]
SAMPLE_INFO = {
    sample["sample_name"]: sample for sample in INPUT_CONTRACT["samples"]
}


def contract_raw_read(sample, mate):
    return SAMPLE_INFO[sample]["reads"][mate]["logical_path"]


def contract_sample_receipt(sample):
    return SAMPLE_INFO[sample]["receipt_path"]


def contract_safe_code(sample):
    return SAMPLE_INFO[sample]["safe_code"]


SPLIT = ["atac", "ct"]
RAW_R2_UNCLIPPED = "processed/{sample}/{split}/{sample}.{split}.R2_5.unclipped.bed.gz"
R2_UNCLIPPED_PHASE = "processed/{sample}/{split}/{sample}.{split}.R2_5.unclipped.phase.tsv.gz"

CUTTAG_FRAGMENT_CONFIG = config.get("cuttag_fragments", {})
if IS_CHARM:
    expected_fragment_contract = {
        "definition": "r2_5prime_unclipped_v2_no_hickit_gate",
        "min_mapq": 20,
        "min_baseq": 20,
        "require_original_r2_qstart_bp": 0,
        "allow_opposite_end_clipping": True,
        "phase_policy": "direct_r2_then_first_phased_r1",
        "endpoint_dedup_eps_bp": 1,
    }
    for key, expected in expected_fragment_contract.items():
        observed = CUTTAG_FRAGMENT_CONFIG.get(key)
        if observed != expected:
            raise ValueError(
                "cuttag_fragments.{} must be {!r}, observed {!r}".format(
                    key, expected, observed
                )
            )

RNA_OUTPUT_TYPE_CHOICES = (
    "r1_all",
    "r1_compatible",
    "r1r2_concordant",
)
RNA_OUTPUT_TYPES = config.get("rna_output_types", [])
if not isinstance(RNA_OUTPUT_TYPES, list) or not RNA_OUTPUT_TYPES:
    raise ValueError(
        "rna_output_types must be a non-empty list chosen from {}".format(
            ", ".join(RNA_OUTPUT_TYPE_CHOICES)
        )
    )
RNA_OUTPUT_TYPES = [str(value).strip() for value in RNA_OUTPUT_TYPES]
if len(RNA_OUTPUT_TYPES) != len(set(RNA_OUTPUT_TYPES)):
    raise ValueError("rna_output_types must not contain duplicate values")
invalid_rna_output_types = sorted(
    set(RNA_OUTPUT_TYPES) - set(RNA_OUTPUT_TYPE_CHOICES)
)
if invalid_rna_output_types:
    raise ValueError(
        "invalid rna_output_types: {}; allowed values are {}".format(
            ", ".join(invalid_rna_output_types),
            ", ".join(RNA_OUTPUT_TYPE_CHOICES),
        )
    )
config["rna_output_types"] = RNA_OUTPUT_TYPES
RNA_FILTERED_OUTPUT_TYPES = [
    value for value in RNA_OUTPUT_TYPES if value != "r1_all"
]
RNA_COMPATIBILITY_ENABLED = bool(RNA_FILTERED_OUTPUT_TYPES)

RNA_COMPATIBILITY_CONFIG = config.get("rna_gene_compatibility", {})
expected_rna_compatibility_contract = {
    "min_mapq": 30,
    "feature_type": "gene",
    "feature_id": "gene_id",
    "allow_overlapping_gene_loci": True,
    "overlap_connected_components": True,
    "r1_strand": 1,
    "r2_strand": 2,
}
for key, expected in expected_rna_compatibility_contract.items():
    observed = RNA_COMPATIBILITY_CONFIG.get(key)
    if observed != expected:
        raise ValueError(
            "rna_gene_compatibility.{} must be {!r}, observed {!r}".format(
                key, expected, observed
            )
        )

CUTADAPT_4_6 = config.get("softwares", {}).get("cutadapt_4_6", "")
if RNA_COMPATIBILITY_ENABLED and not CUTADAPT_4_6:
    raise ValueError(
        "softwares.cutadapt_4_6 is required for filtered RNA output types"
    )

METADATA_JUPYTER = config.get("softwares", {}).get("metadata_jupyter", "")
METADATA_KERNEL = config.get("metadata", {}).get("kernel_name", "")
if not METADATA_JUPYTER or not METADATA_KERNEL:
    raise ValueError(
        "softwares.metadata_jupyter and metadata.kernel_name are required"
    )

PRIMARY_RNA_OUTPUT_TYPE = str(
    config.get("rna_primary_output_type", "")
).strip()
if PRIMARY_RNA_OUTPUT_TYPE not in RNA_OUTPUT_TYPES:
    raise ValueError(
        "rna_primary_output_type must be one of the selected rna_output_types; "
        "observed {!r}".format(PRIMARY_RNA_OUTPUT_TYPE)
    )
config["rna_primary_output_type"] = PRIMARY_RNA_OUTPUT_TYPE
PRIMARY_RNA_COUNT_BAM = (
    "processed/RNA_all/feature_{}_gene_total/samsort.bam".format(
        PRIMARY_RNA_OUTPUT_TYPE
    )
)

#############TARGET CONTRACTS###############

RNA_GENOMES = ["total", "genome1", "genome2"] if config["if_RNA_snp_split"] else ["total"]
STRUCTURE_RESOLUTIONS = ["20k", "50k", "200k", "1m"]
STRUCTURE_REPLICATES = list(range(5))


def _unique_paths(paths):
    seen = set()
    unique = []
    for path in paths:
        if path not in seen:
            seen.add(path)
            unique.append(path)
    return unique


RNA_MATRIX_TARGET_OUTPUTS = _unique_paths(
    expand(
        "result/RNA_Res/{rna_output_type}/counts.{type}.{genome}.safe.tsv",
        rna_output_type=RNA_OUTPUT_TYPES,
        type=["gene", "exon"],
        genome=RNA_GENOMES,
    )
    + expand(
        "result/RNA_Res/{rna_output_type}/counts.{type}.{genome}.tsv",
        rna_output_type=RNA_OUTPUT_TYPES,
        type=["gene", "exon"],
        genome=RNA_GENOMES,
    )
    + expand(
        "result/RNA_Res/{rna_output_type}/counts.{type}.{genome}.format.tsv",
        rna_output_type=RNA_OUTPUT_TYPES,
        type=["gene", "exon"],
        genome=RNA_GENOMES,
    )
    + expand(
        "result/RNA_Res/{rna_output_type}/counts.{type}.{genome}.cell_contract.tsv",
        rna_output_type=RNA_OUTPUT_TYPES,
        type=["gene", "exon"],
        genome=RNA_GENOMES,
    )
)
RNA_MODE_QC_OUTPUTS = [
    "qc/stat/rna.output_modes.per_cell.tsv",
    "qc/stat/rna.output_modes.summary.tsv",
]
RNA_COMPATIBILITY_QC_OUTPUTS = (
    [
        "qc/stat/rna.r2_polyT.cutadapt.json",
        "qc/stat/rna.r2_validator.star.Log.final.out",
        "qc/stat/rna.gene_compatibility.per_cell.tsv",
    ]
    if RNA_COMPATIBILITY_ENABLED
    else []
)
RNA_QC_OUTPUTS = RNA_MODE_QC_OUTPUTS + RNA_COMPATIBILITY_QC_OUTPUTS
RNA_TARGET_OUTPUTS = _unique_paths(
    expand("processed/{sample}/umi/umi.{sample}.rna.R2.fq.gz", sample=SAMPLES)
    + RNA_MATRIX_TARGET_OUTPUTS
    + RNA_QC_OUTPUTS
)

CONTACTS2D_TARGET_OUTPUTS = _unique_paths(
    expand("processed/{sample}/2d_info/contacts.pairs.gz", sample=SAMPLES)
    + expand("result/cleaned_pairs/c12/{sample}.pairs.gz", sample=SAMPLES)
    + expand("result/cleaned_pairs/c123/{sample}.pairs.gz", sample=SAMPLES)
)

FRAGMENT_TARGET_OUTPUTS = _unique_paths(
    expand(
        "processed/{sample}/{split}/{sample}.{split}.R1.fq.gz",
        sample=SAMPLES if IS_CHARM else [],
        split=SPLIT if IS_CHARM else [],
    )
    + expand(
        RAW_R2_UNCLIPPED,
        sample=SAMPLES if IS_CHARM else [],
        split=SPLIT if IS_CHARM else [],
    )
    + expand(
        R2_UNCLIPPED_PHASE,
        sample=SAMPLES if IS_CHARM else [],
        split=SPLIT if IS_CHARM else [],
    )
    + expand(
        "result/fragments/{split}.fragments.bgz",
        split=SPLIT if IS_CHARM else [],
    )
    + expand(
        "result/fragments/{split}.fragments.bgz.tbi",
        split=SPLIT if IS_CHARM else [],
    )
)

STRUCTURE3D_TARGET_OUTPUTS = _unique_paths(
    expand(
        "result/3d_info/{sample}/{sample}.{res}.align.rms.info",
        sample=SAMPLES if config["if_structure"] else [],
        res=STRUCTURE_RESOLUTIONS if config["if_structure"] else [],
    )
    + expand(
        "result/3d_info/{sample}/clean.{res}.{rep}.3dg.gz",
        sample=SAMPLES if config["if_structure"] else [],
        res=STRUCTURE_RESOLUTIONS if config["if_structure"] else [],
        rep=STRUCTURE_REPLICATES if config["if_structure"] else [],
    )
)

QC_COMMON_STAT_OUTPUTS = [
    "qc/stat/raw.pairs.stat",
    "qc/stat/raw.fq.stat",
    "qc/stat/rna.fq.stat",
    "qc/stat/dna.fq.stat",
    "qc/stat/pairs.dedup.stat",
    "qc/stat/pairs.c1.stat",
    "qc/stat/pairs.c12.stat",
    "qc/stat/pairs.c123.stat",
    "qc/stat/inter.pairs.c123.stat",
    "qc/stat/yperx.stat",
    "qc/stat/rna.reads_per_cell.stat",
    "qc/stat/rna.dna_contam.stat",
]
QC_CHARM_STAT_OUTPUTS = (
    [
        "qc/stat/atac.read.stat",
        "qc/stat/ct.read.stat",
        "qc/stat/atac.dedup_rate.stat",
        "qc/stat/ct.dedup_rate.stat",
    ]
    if IS_CHARM
    else []
)
QC_STRUCTURE_STAT_OUTPUTS = (
    ["qc/stat/rmsd.info"]
    if config["if_structure"]
    else []
)
QC_STAT_OUTPUTS = (
    QC_COMMON_STAT_OUTPUTS + QC_CHARM_STAT_OUTPUTS + QC_STRUCTURE_STAT_OUTPUTS
)
NOTEBOOK_OUTPUTS = ["qc/stat.ipynb", "qc/metadata_qc.R"]
METADATA_OUTPUTS = ["qc/stat.executed.ipynb", "qc/metadata_raw.tsv"]
AUDIT_OUTPUTS = ["qc/COMPLETE_RUN_AUDIT.tsv"]
PROVENANCE_OUTPUTS = [
    "qc/provenance/effective_config.json",
    "qc/provenance/source_files.sha256.tsv",
]
TARGET_CONTRACT_OUTPUTS = ["qc/target_outputs.tsv"]
QC_BASE_OUTPUTS = _unique_paths(
    QC_STAT_OUTPUTS
    + RNA_QC_OUTPUTS
    + NOTEBOOK_OUTPUTS
    + PROVENANCE_OUTPUTS
    + TARGET_CONTRACT_OUTPUTS
)
CORE_TARGET_OUTPUTS = _unique_paths(
    RNA_TARGET_OUTPUTS + CONTACTS2D_TARGET_OUTPUTS + FRAGMENT_TARGET_OUTPUTS
)
DELIVERY_TARGET_OUTPUTS = _unique_paths(
    CORE_TARGET_OUTPUTS
    + STRUCTURE3D_TARGET_OUTPUTS
    + QC_BASE_OUTPUTS
    + METADATA_OUTPUTS
)
QC_TARGET_OUTPUTS = _unique_paths(QC_BASE_OUTPUTS + METADATA_OUTPUTS + AUDIT_OUTPUTS)
ALL_TARGET_OUTPUTS = _unique_paths(DELIVERY_TARGET_OUTPUTS + AUDIT_OUTPUTS)

TARGET_GROUPS = [
    ("rna", True, RNA_TARGET_OUTPUTS),
    ("contacts2d", True, CONTACTS2D_TARGET_OUTPUTS),
    ("fragments", IS_CHARM, FRAGMENT_TARGET_OUTPUTS),
    ("structure3d", bool(config["if_structure"]), STRUCTURE3D_TARGET_OUTPUTS),
    ("provenance", True, PROVENANCE_OUTPUTS),
    ("qc", True, QC_TARGET_OUTPUTS),
    ("all", True, ALL_TARGET_OUTPUTS),
]


def _render_target_contract():
    lines = ["target\tenabled\toutput"]
    for target_name, enabled, paths in TARGET_GROUPS:
        if not paths:
            lines.append("{}\t{}\tNA".format(target_name, int(enabled)))
            continue
        lines.extend(
            "{}\t{}\t{}".format(target_name, int(enabled), path)
            for path in paths
        )
    return "\n".join(lines) + "\n"


EFFECTIVE_CONFIG_TEXT = json.dumps(dict(config), indent=2, sort_keys=True) + "\n"
EFFECTIVE_CONFIG_SHA256 = hashlib.sha256(
    EFFECTIVE_CONFIG_TEXT.encode("utf-8")
).hexdigest()
TARGET_CONTRACT_TEXT = _render_target_contract()
TARGET_CONTRACT_SHA256 = hashlib.sha256(
    TARGET_CONTRACT_TEXT.encode("utf-8")
).hexdigest()

WORKFLOW_SOURCE_INVENTORY = build_source_inventory(
    PIPELINE_DIR, CONFIG_SOURCE_FILES, CHARMTOOLS_DIR
)
WORKFLOW_SOURCE_FILES = [path for _, path in WORKFLOW_SOURCE_INVENTORY]

QC_UPSTREAM_OUTPUTS = _unique_paths(
    CORE_TARGET_OUTPUTS + STRUCTURE3D_TARGET_OUTPUTS
)
QC_RAW_FASTQ_INPUTS = [
    contract_raw_read(sample, mate)
    for sample in SAMPLES
    for mate in ("r1", "r2")
]
METADATA_UPSTREAM_OUTPUTS = _unique_paths(
    QC_STAT_OUTPUTS
    + RNA_QC_OUTPUTS
    + RNA_MATRIX_TARGET_OUTPUTS
    + NOTEBOOK_OUTPUTS
    + PROVENANCE_OUTPUTS
    + expand(
        "result/fragments/{split}.fragments.bgz",
        split=SPLIT if IS_CHARM else [],
    )
    + expand(
        "result/fragments/{split}.fragments.bgz.tbi",
        split=SPLIT if IS_CHARM else [],
    )
)


#############DECLARATIVE TARGETS###############

rule all:
    input:
        ALL_TARGET_OUTPUTS


rule rna:
    input:
        RNA_TARGET_OUTPUTS


rule contacts2d:
    input:
        CONTACTS2D_TARGET_OUTPUTS


rule fragments:
    input:
        FRAGMENT_TARGET_OUTPUTS


if config["if_structure"]:
    rule structure3d:
        input:
            STRUCTURE3D_TARGET_OUTPUTS


rule qc:
    input:
        QC_TARGET_OUTPUTS


rule provenance:
    input:
        PROVENANCE_OUTPUTS


rule generate_statistics:
    input:
        input_contract = INPUT_CONTRACT_PATH,
        rna_count_bam = PRIMARY_RNA_COUNT_BAM,
        raw_fastqs = QC_RAW_FASTQ_INPUTS,
        upstream = QC_UPSTREAM_OUTPUTS,
        script = os.path.join(SCRIPTS_DIR, "generateStat.sh"),
        contract_script = os.path.join(SCRIPTS_DIR, "generate_stat_contract.py"),
    output:
        QC_STAT_OUTPUTS
    threads: config["resources"]["generateStat_cpu_threads"]
    params:
        structure_enabled = str(bool(config["if_structure"])),
        experiment_type = EXPERIMENT_TYPE,
        structure_resolutions = ",".join(STRUCTURE_RESOLUTIONS),
    shell:
        "CHARM_STAT_THREADS={threads} bash {input.script} {WORK_DIR} "
        "{input.input_contract} {params.structure_enabled} "
        "{params.experiment_type} {params.structure_resolutions} "
        "{input.rna_count_bam}"


rule stage_qc_notebook:
    input:
        notebook = os.path.join(SCRIPTS_DIR, "stat.ipynb"),
        helper = os.path.join(SCRIPTS_DIR, "metadata_qc.R"),
    output:
        notebook = NOTEBOOK_OUTPUTS[0],
        helper = NOTEBOOK_OUTPUTS[1],
    shell:
        "mkdir -p qc && cp {input.notebook} {output.notebook} && "
        "cp {input.helper} {output.helper}"


rule generate_metadata:
    input:
        notebook = NOTEBOOK_OUTPUTS[0],
        jupyter = METADATA_JUPYTER,
        upstream = METADATA_UPSTREAM_OUTPUTS,
        reference = (
            [config["refs"][config["ref_genome"]]["bwa_mem2_index"] + ".fai",
             config["refs"][config["ref_genome"]]["annotations"]]
            if IS_CHARM and config["ref_genome"] != "GRCm38" else []
        ),
    output:
        executed_notebook = METADATA_OUTPUTS[0],
        metadata = METADATA_OUTPUTS[1],
    log:
        main = "qc/logs/metadata.log",
    params:
        kernel = METADATA_KERNEL,
    shell:
        r"""
        set -euo pipefail
        mkdir -p qc/logs tmp/ipython
        export IPYTHONDIR="{WORK_DIR}/tmp/ipython"
        (
            cd qc
            export PATH="$(dirname "{input.jupyter}"):$PATH"
            {input.jupyter} nbconvert --to notebook --execute stat.ipynb \
                --output stat.executed.ipynb \
                --ExecutePreprocessor.timeout=-1 \
                --ExecutePreprocessor.kernel_name={params.kernel}
            test -s metadata_raw.tsv
        ) > {log.main} 2>&1
        """


rule audit_complete_run:
    input:
        deliverables = DELIVERY_TARGET_OUTPUTS,
        script = os.path.join(SCRIPTS_DIR, "audit_complete_run.py"),
        contract_script = os.path.join(SCRIPTS_DIR, "generate_stat_contract.py"),
    output:
        receipt = AUDIT_OUTPUTS[0],
    log:
        main = "qc/logs/audit.log",
    shell:
        r"""
        set -euo pipefail
        {{
            python {input.script} --work-dir {WORK_DIR} --pipeline-dir {PIPELINE_DIR}

            find processed -type f -name '*.bam' -print0 > tmp/audit_bam_files
            bam_count=$(tr -cd '\0' < tmp/audit_bam_files | wc -c)
            xargs -0 -r samtools quickcheck < tmp/audit_bam_files
            printf 'bam_quickcheck\tPASS\t%s BAM files\n' "$bam_count" >> {output.receipt}

            if [ "{EXPERIMENT_TYPE}" = charm ]; then
                tabix -l result/fragments/atac.fragments.bgz >/dev/null
                tabix -l result/fragments/ct.fragments.bgz >/dev/null
                printf 'fragment_indexes\tPASS\tatac and ct\n' >> {output.receipt}
            fi

            if find -L Rawdata processed result qc -type l -print -quit | grep -q .; then
                printf 'broken symlink detected\n' >&2
                exit 2
            fi
            printf 'symlink_integrity\tPASS\tno broken links\n' >> {output.receipt}
        }} > {log.main} 2>&1
        """


rule pipeline_provenance:
    input:
        sources = WORKFLOW_SOURCE_FILES,
    output:
        effective_config = PROVENANCE_OUTPUTS[0],
        source_hashes = PROVENANCE_OUTPUTS[1],
    params:
        effective_config_sha256 = EFFECTIVE_CONFIG_SHA256,
    run:
        os.makedirs(os.path.dirname(str(output.effective_config)), exist_ok=True)
        with open(str(output.effective_config), "w") as handle:
            handle.write(EFFECTIVE_CONFIG_TEXT)
        write_source_hash_inventory(
            WORKFLOW_SOURCE_INVENTORY, PIPELINE_DIR, str(output.source_hashes)
        )


rule target_contract:
    input:
        snakefile = os.path.join(PIPELINE_DIR, "CHARM.smk"),
        config_sources = CONFIG_SOURCE_FILES,
    output:
        contract = TARGET_CONTRACT_OUTPUTS[0],
    params:
        contract_sha256 = TARGET_CONTRACT_SHA256,
    run:
        os.makedirs(os.path.dirname(str(output.contract)), exist_ok=True)
        with open(str(output.contract), "w") as handle:
            handle.write(TARGET_CONTRACT_TEXT)

    
############END_rule_all############


include: os.path.join(PIPELINE_DIR, "rules/CHARM_split.rules")
include: os.path.join(PIPELINE_DIR, "rules/CHARM_cuttag.rules")
include: os.path.join(PIPELINE_DIR, "rules/scHiC_2dprocess.rules")
if config["if_structure"]:
    include: os.path.join(PIPELINE_DIR, "rules/scHiC_3dprocess.rules")
include: os.path.join(PIPELINE_DIR, "rules/CHARM_RNA.rules")

from snakemake_compat import track_rule_changes
track_rule_changes(workflow)

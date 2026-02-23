####################################
#       CHARM_pipline              #
#@author Z Liu                     #
#@Ver 0.3.0                        #
#@date 2023/8/7                    #
####################################

#############CONFIG#################

import os
from glob import glob

PIPELINE_DIR = workflow.basedir
SCRIPTS_DIR = os.path.join(PIPELINE_DIR, "CHARM_scripts")
RUN_ROOT = os.path.dirname(PIPELINE_DIR)

CONFIG_PATH = os.path.join(PIPELINE_DIR, "config.yaml")
configfile: CONFIG_PATH

CHARMTOOLS_DIR = config.get("softwares", {}).get("CHARMtools", os.path.join(RUN_ROOT, "CHARMtools"))
if not os.path.isabs(CHARMTOOLS_DIR):
    CHARMTOOLS_DIR = os.path.normpath(os.path.join(PIPELINE_DIR, CHARMTOOLS_DIR))
CHARMTOOLS_PYTHONPATH = os.path.dirname(CHARMTOOLS_DIR)

WORK_DIR = config.get("work_dir", "..")
if isinstance(WORK_DIR, str):
    WORK_DIR = WORK_DIR.replace("\r", "").strip()
if not os.path.isabs(WORK_DIR):
    WORK_DIR = os.path.normpath(os.path.join(PIPELINE_DIR, WORK_DIR))
workdir: WORK_DIR

# input
def _discover_samples(raw_dir):
    samples = []
    for entry in os.listdir(raw_dir):
        entry_path = os.path.join(raw_dir, entry)
        if not os.path.isdir(entry_path):
            continue
        candidates = [
            (f"{entry}_1.fq.gz", f"{entry}_2.fq.gz"),
            (f"{entry}_R1.fq.gz", f"{entry}_R2.fq.gz"),
        ]
        for r1_name, r2_name in candidates:
            r1 = os.path.join(entry_path, r1_name)
            r2 = os.path.join(entry_path, r2_name)
            if os.path.exists(r1) and os.path.exists(r2):
                # treat empty fastq as missing
                if os.path.getsize(r1) > 0 and os.path.getsize(r2) > 0:
                    samples.append(entry)
                break
    return sorted(samples)

SAMPLES = _discover_samples(os.path.join(WORK_DIR, "Rawdata"))
SPLIT = ["atac", "ct"]

#############RULE_ALL###############
"""
decide what you need for your down stream analysis.
"""
RELAXED_RULE_ALL = bool(config.get("relaxed_rule_all", False))

def _rule_all_inputs():
    if not RELAXED_RULE_ALL:
        return [
            #preliminary split
            expand("processed/{sample}/umi/umi.{sample}.rna.R2.fq.gz",sample=SAMPLES),
            #RNA part
            expand("result/RNA_Res/counts.{type}.{genome}.tsv",type=["gene","exon"],genome=["total","genome1","genome2"] if config["if_RNA_snp_split"] else ["total"]),
            expand("result/RNA_Res/counts.{type}.{genome}.format.tsv",type=["gene","exon"],genome=["total","genome1","genome2"] if config["if_RNA_snp_split"] else ["total"]),
            #Hi-C part pairs info
            expand("processed/{sample}/2d_info/contacts.pairs.gz",sample=SAMPLES),
            expand("result/cleaned_pairs/c12/{sample}.pairs.gz",sample=SAMPLES),
            #Hi-C part 3d info
            expand("result/3d_info/{sample}/{sample}.{res}.align.rms.info",sample=SAMPLES if config["if_structure"] else [],res=["20k","50k","200k","1m"] if config["if_structure"] else []),
            expand("result/3d_info/{sample}/clean.{res}.{rep}.3dg.gz", sample=SAMPLES if config["if_structure"] else [],
                res=["20k","50k","200k","1m"] if config["if_structure"] else [],
                rep=list(range(5)) if config["if_structure"] else []),

            #cuttag part
            expand("processed/{sample}/{split}/{sample}.{split}.R1.fq.gz", sample=SAMPLES if config["if_charm"] else [],split=SPLIT if config ["if_charm"] else []),
            expand("processed/{sample}/{split}/{sample}.{split}.R2_5.bed.gz", sample=SAMPLES if config["if_charm"] else [],split=SPLIT if config ["if_charm"] else []),
            expand("result/fragments/{split}.fragments.bgz",split=SPLIT if config ["if_charm"] else [])
        ]

    # Relaxed mode: only require outputs that already exist.
    inputs = []
    # preliminary split (RNA UMI)
    inputs += glob(os.path.join(WORK_DIR, "processed/*/umi/umi.*.rna.R2.fq.gz"))
    # RNA part
    inputs += glob(os.path.join(WORK_DIR, "result/RNA_Res/counts.*.tsv"))
    inputs += glob(os.path.join(WORK_DIR, "result/RNA_Res/counts.*.format.tsv"))
    # Hi-C pairs
    inputs += glob(os.path.join(WORK_DIR, "processed/*/2d_info/contacts.pairs.gz"))
    inputs += glob(os.path.join(WORK_DIR, "result/cleaned_pairs/c12/*.pairs.gz"))
    # Hi-C 3D
    if config["if_structure"]:
        inputs += glob(os.path.join(WORK_DIR, "result/3d_info/*/*.align.rms.info"))
        inputs += glob(os.path.join(WORK_DIR, "result/3d_info/*/clean.*.3dg.gz"))
    # cuttag part
    if config["if_charm"]:
        for split in SPLIT:
            inputs += glob(os.path.join(WORK_DIR, f"processed/*/{split}/*.{split}.R1.fq.gz"))
            inputs += glob(os.path.join(WORK_DIR, f"processed/*/{split}/*.{split}.R2_5.bed.gz"))
            frag = os.path.join(WORK_DIR, f"result/fragments/{split}.fragments.bgz")
            if os.path.exists(frag):
                inputs.append(frag)
    return inputs

rule all:
    input:
        _rule_all_inputs()

    threads: config["resources"]["generateStat_cpu_threads"]
    shell:"""
        {SCRIPTS_DIR}/generateStat.sh
        mkdir -p {WORK_DIR}/analysis
        cp {SCRIPTS_DIR}/stat.ipynb {WORK_DIR}/analysis/
        echo "done!"
    """

    
############END_rule_all############


include: os.path.join(PIPELINE_DIR, "rules/CHARM_split.rules")
include: os.path.join(PIPELINE_DIR, "rules/CHARM_cuttag.rules")
include: os.path.join(PIPELINE_DIR, "rules/scHiC_2dprocess.rules")
include: os.path.join(PIPELINE_DIR, "rules/scHiC_3dprocess.rules")
include: os.path.join(PIPELINE_DIR, "rules/CHARM_RNA.rules")

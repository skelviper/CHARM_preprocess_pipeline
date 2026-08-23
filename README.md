# CHARM Preprocess Pipeline

## Purpose

This repository contains the Snakemake preprocessing pipeline for single-cell CHARM and HiRES data. Starting from paired-end FASTQ files for each cell, it produces RNA expression matrices, DNA 2D contacts, CHARM ATAC/CUT&Tag fragments, and consolidated QC metadata.

The pipeline entry point is `runCHARM.sh`. The run directory is selected by `work_dir` in `config.yaml`.

## Input Layout

The `work_dir` directory must contain `Rawdata/`. Each cell has its own directory, and the directory name is used as `<sample>`.

Two paired-end FASTQ naming schemes are supported: `<sample>_1.fq.gz` + `<sample>_2.fq.gz`, or `<sample>_R1.fq.gz` + `<sample>_R2.fq.gz`.

```text
<work_dir>/
└── Rawdata/
    ├── <sample_A>/
    │   ├── <sample_A>_1.fq.gz
    │   └── <sample_A>_2.fq.gz
    └── <sample_B>/
        ├── <sample_B>_R1.fq.gz
        └── <sample_B>_R2.fq.gz
```

Each cell must contain exactly one complete R1/R2 pair. The two naming schemes cannot be mixed within one cell directory, and neither mate may be missing.

## Output Layout

All outputs are written under `work_dir`:

```text
<work_dir>/
├── Rawdata/                         # Input FASTQ files
├── processed/
│   └── <sample>/                    # Per-cell FASTQ, BAM, and intermediate files
├── result/
│   ├── RNA_Res/
│   │   └── <configured_output>/     # Gene and exon RNA matrices
│   ├── cleaned_pairs/
│   │   ├── c1/                      # Pairs after promiscuous-leg removal
│   │   ├── c12/                     # Pairs after isolated-contact removal
│   │   └── c123/                    # Pairs after splicing-contact removal
│   ├── fragments/                   # CHARM ATAC/CUT&Tag fragments and indexes
│   └── 3d_info/                     # Optional 3D structures
├── qc/
│   ├── metadata_raw.tsv             # Main per-cell metadata table
│   ├── stat/                        # Raw summary tables
│   ├── logs/                        # Workflow and rule logs
│   ├── input_contract/              # Input list used by the run
│   ├── provenance/                  # Effective configuration and source record
│   ├── stat.executed.ipynb          # Executed statistics notebook
│   ├── target_outputs.tsv           # Outputs required by the active configuration
│   └── COMPLETE_RUN_AUDIT.tsv       # Final complete-run audit
└── tmp/                              # Runtime temporary files; removed on launcher exit
```

`processed/` contains per-cell working files and is mainly useful for troubleshooting. Routine analysis generally uses the matrices, pairs, and fragments under `result/`, together with `qc/metadata_raw.tsv` and the summary tables under `qc/stat/`. Optional directories depend on `experiment_type` and `if_structure`.

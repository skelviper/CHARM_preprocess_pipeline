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
└── tmp/                              # Runtime temporary files; each launcher cleans its own directory
```

`processed/` contains per-cell working files and is mainly useful for troubleshooting. Routine analysis generally uses the matrices, pairs, and fragments under `result/`, together with `qc/metadata_raw.tsv` and the summary tables under `qc/stat/`. Optional directories depend on `experiment_type` and `if_structure`.

## Resuming and QC

The pinned Snakemake 5.20 runtime now schedules jobs when their recorded rule
parameters, input paths, or rule code change. This also applies to direct
`snakemake -s CHARM.smk` invocations. Keep the run's `.snakemake/metadata/`
directory so these comparisons remain available. A dry-run reports the work
without deleting existing temporary files.

ATAC/CT read estimates exclude secondary and supplementary alignments. The
existing metadata conversion still reports gigabases assuming 300 bp per read
pair. `nCount_atac` and `nCount_ct` remain 5 kb bin counts, not unique fragment
counts.

CHARM TSS QC uses the original mm10/EnsDb v79 annotation for `GRCm38`. Other
references, including `GRCh37d5`, use the configured `refs` FASTA `.fai` and GTF
transcripts with `gene_id`, `gene_name`, and `gene_biotype` or `gene_type`.
Chromosome names must match the FASTA. TSS scores retain Signac's original
`fast = FALSE` pileup, center/flank windows, and background normalization;
position-enrichment matrices are not stored in a Seurat object.

An RNA matrix with no counted UMIs retains every frozen cell column and zero
feature rows. Fragment QC runs independently of RNA. Unmeasurable TSS scores
are `NA` with `TSS.status.atac` or `TSS.status.ct` set to `insufficient_data`.
Enabled 3D outputs with no common particles contain
`#status\tinsufficient_data`; their `rmsd_*` metadata values remain `NA`.

#!/bin/bash

#usage: ./runCHARM.sh

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${PIPELINE_DIR}/config.yaml"
WORK_DIR="$(awk -F': *' '/^work_dir:/ {print $2; exit}' "${CONFIG_FILE}")"
if [ -z "${WORK_DIR}" ]; then
  WORK_DIR=".."
fi
if [[ "${WORK_DIR}" != /* ]]; then
  WORK_DIR="${PIPELINE_DIR}/${WORK_DIR}"
fi
WORK_DIR="$(cd "${WORK_DIR}" && pwd)"

mkdir -p "${WORK_DIR}/slurm_log"
cd "${WORK_DIR}"
snakemake --use-conda --cluster 'sbatch --qos=high -w node03 --output=slurm_log/slurm-%j.out --cpus-per-task={threads} -t 7-00:00 -J CHARM!' --jobs 1024 --rerun-incomplete -s "${PIPELINE_DIR}/CHARM.smk" --keep-going

mkdir -p "${WORK_DIR}/analysis"
cp "${PIPELINE_DIR}/stat.ipynb" "${WORK_DIR}/analysis/stat.ipynb"

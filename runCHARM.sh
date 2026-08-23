#!/usr/bin/env bash

set -u

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CHARM_ENV_NAME="charm"
CHARM_PYTHON_BIN="${CHARM_PYTHON_BIN:-python3}"

usage() {
    cat <<'EOF'
Usage: ./runCHARM.sh [--accept-input-change] [SNAKEMAKE_ARGUMENT ...]

The first launch freezes automatically discovered Rawdata inputs. Later launches
validate that frozen contract without rewriting it. Use --accept-input-change
only for an intentional cohort/path replacement; the prior generated contract
is archived before the replacement is published.

runCHARM.sh resolves config.yaml plus Snakemake --configfile and --config
work_dir=... overrides before touching run state. It owns --snakefile and
--directory and rejects --profile. It runs the entire workflow in the existing
charm Conda environment, uses work_dir/tmp, and removes that temporary tree
when the launcher exits.
Mapping uses one BWA job per cell. Jobs are
submitted directly through Snakemake's --cluster interface, matching the
frozen Part 1 launch model. Real Slurm and complete 200-cell execution gates
both remain NOT_RUN.
EOF
}

for argument in "$@"; do
    if [ "${argument}" = "-h" ] || [ "${argument}" = "--help" ]; then
        usage
        exit 0
    fi
done

if [ "${CONDA_DEFAULT_ENV:-}" != "${CHARM_ENV_NAME}" ]; then
    CHARM_CONDA_BIN="${CHARM_CONDA_BIN:-conda}"
    if ! command -v "${CHARM_CONDA_BIN}" >/dev/null 2>&1; then
        echo "ERROR: conda is required to enter the existing ${CHARM_ENV_NAME} environment" >&2
        exit 2
    fi
    exec "${CHARM_CONDA_BIN}" run --no-capture-output -n "${CHARM_ENV_NAME}" "$0" "$@"
fi

accept_input_change=0
snakemake_args=()
while [ "$#" -gt 0 ]; do
    case "$1" in
        --accept-input-change)
            if [ "${accept_input_change}" -eq 1 ]; then
                echo "ERROR: --accept-input-change was provided more than once" >&2
                exit 2
            fi
            accept_input_change=1
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            snakemake_args+=("$1")
            ;;
    esac
    shift
done

INVOCATION_DIR="${PWD}"
SNAKEMAKE_BIN="${CHARM_SNAKEMAKE_BIN:-snakemake}"
for required_command in "${CHARM_PYTHON_BIN}" "${SNAKEMAKE_BIN}"; do
    if ! command -v "${required_command}" >/dev/null 2>&1; then
        echo "ERROR: required launcher executable is unavailable: ${required_command}" >&2
        exit 2
    fi
done

resolution_root="${TMPDIR:-${INVOCATION_DIR}}"
mkdir -p "${resolution_root}"
resolution_file="$(mktemp "${resolution_root%/}/.charm-run-config.XXXXXX")" || exit 2
if ! "${CHARM_PYTHON_BIN}" "${PIPELINE_DIR}/CHARM_scripts/resolve_run_config.py" \
    --pipeline-dir "${PIPELINE_DIR}" \
    --invocation-dir "${INVOCATION_DIR}" \
    -- "${snakemake_args[@]}" > "${resolution_file}"; then
    rm -f "${resolution_file}"
    echo "Effective run configuration could not be resolved; no run state was changed." >&2
    exit 2
fi
mapfile -d '' -t resolved_run < "${resolution_file}"
rm -f "${resolution_file}"
if [ "${#resolved_run[@]}" -lt 1 ]; then
    echo "ERROR: config resolver returned no effective work_dir" >&2
    exit 2
fi
WORK_DIR="${resolved_run[0]}"
snakemake_args=("${resolved_run[@]:1}")
export CHARM_EFFECTIVE_WORK_DIR="${WORK_DIR}"

WORKSPACE_TMP="${WORK_DIR}/tmp"
if [ -z "${TMPDIR:-}" ]; then
    TMPDIR="${WORKSPACE_TMP}"
elif [[ "${TMPDIR}" != /* ]]; then
    TMPDIR="${INVOCATION_DIR}/${TMPDIR}"
fi
mkdir -p "${TMPDIR}"
TMPDIR="$(cd "${TMPDIR}" && pwd)"
export TMPDIR

cleanup_workspace_tmp() {
    rm -rf -- "${WORKSPACE_TMP}"
}
trap cleanup_workspace_tmp EXIT

contract_state="${WORK_DIR}/qc/input_contract"
contract_action="create"
if [ "${accept_input_change}" -eq 1 ]; then
    contract_action="replace"
elif [ -e "${contract_state}" ] || [ -L "${contract_state}" ]; then
    contract_action="validate"
fi
if ! "${CHARM_PYTHON_BIN}" "${PIPELINE_DIR}/CHARM_scripts/input_contract.py" \
    "${contract_action}" --work-dir "${WORK_DIR}"; then
    echo "Rawdata input-contract ${contract_action} failed; workflow was not submitted." >&2
    exit 2
fi

mkdir -p "${WORK_DIR}/qc/logs"
cd "${WORK_DIR}"
"${SNAKEMAKE_BIN}" \
    --cluster 'sbatch --qos=high -w node03 --output=qc/logs/slurm-%j.out --cpus-per-task={threads} -t 7-00:00 -J CHARM!' \
    --jobs 1024 \
    --resources star_slots=1 count_slots=1 \
    --rerun-incomplete \
    -s "${PIPELINE_DIR}/CHARM.smk" \
    "${snakemake_args[@]}"

#!/bin/bash
# Shared by the UQ job scripts (sourced, not submitted).
#
#  * activates the gwas2m-lncrna environment
#  * picks the configuration: $GWAS2M_LNCRNA_CONFIG, else ./gwas2m-lncrna.yaml in
#    the submit folder, else the UQ sample next to this file
#  * run_with_resubmit: stops the command before the 24 h UQ limit and resubmits
#    the same job; finished work is skipped, so the next job continues where
#    this one stopped.

set -eo pipefail

# How long one job may run before it stops itself (headroom below 24:00:00).
GWAS2M_TIME_LIMIT="${GWAS2M_TIME_LIMIT:-23h}"
# Maximum number of automatic resubmissions per run.
GWAS2M_MAX_RESUBMIT="${GWAS2M_MAX_RESUBMIT:-10}"
GWAS2M_ATTEMPT="${GWAS2M_ATTEMPT:-1}"

cd "${SLURM_SUBMIT_DIR:-$PWD}"

# Activate before `set -u`: conda's activation scripts use unset variables.
# Skipped when the environment is already active (sbatch passes your shell's
# environment to the job, so submitting from an activated shell is enough).
if ! command -v gwas2m-lncrna >/dev/null 2>&1; then
    if command -v conda >/dev/null 2>&1; then
        eval "$(conda shell.bash hook)"
        conda activate gwas2m-lncrna
    elif command -v mamba >/dev/null 2>&1; then
        eval "$(mamba shell hook --shell bash)"
        mamba activate gwas2m-lncrna
    else
        echo "Neither conda nor mamba found; load your conda module or activate gwas2m-lncrna before sbatch." >&2
        exit 1
    fi
fi
if ! command -v gwas2m-lncrna >/dev/null 2>&1; then
    echo "gwas2m-lncrna is not available after activation; see 'Install troubleshooting' in the README." >&2
    exit 1
fi
set -u

if [[ -n "${GWAS2M_LNCRNA_CONFIG:-}" ]]; then
    CONFIG="$GWAS2M_LNCRNA_CONFIG"
elif [[ -f gwas2m-lncrna.yaml ]]; then
    CONFIG="$PWD/gwas2m-lncrna.yaml"
else
    CONFIG="$UQ_DIR/gwas2m-lncrna.uq.yaml"
fi
export GWAS2M_LNCRNA_CONFIG="$CONFIG"

echo "Job          : ${SLURM_JOB_ID:-local} (attempt $GWAS2M_ATTEMPT of $((GWAS2M_MAX_RESUBMIT + 1)))"
echo "Host         : $(hostname)"
echo "Working dir  : $PWD"
echo "Configuration: $CONFIG"

run_with_resubmit() {
    local status=0
    timeout --signal=TERM --kill-after=10m "$GWAS2M_TIME_LIMIT" "$@" || status=$?
    if [[ $status -eq 124 || $status -eq 137 ]]; then
        if [[ -z "${SLURM_JOB_ID:-}" ]]; then
            echo "Time limit reached (not running under SLURM, so not resubmitting)." >&2
            return 1
        fi
        if (( GWAS2M_ATTEMPT > GWAS2M_MAX_RESUBMIT )); then
            echo "Time limit reached and GWAS2M_MAX_RESUBMIT=$GWAS2M_MAX_RESUBMIT used up; resubmit manually." >&2
            return 1
        fi
        local array=()
        [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] && array=(--array="$SLURM_ARRAY_TASK_ID")
        echo "Time limit reached; resubmitting $SCRIPT_PATH to continue (finished work is kept)."
        sbatch ${array[@]+"${array[@]}"} --export=ALL,GWAS2M_ATTEMPT=$((GWAS2M_ATTEMPT + 1)) "$SCRIPT_PATH"
        return 0
    fi
    return $status
}

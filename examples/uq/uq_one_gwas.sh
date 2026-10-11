#!/bin/bash
# UQ HPC: steps 2-7 for ONE GWAS study; submit as an array, one task per study.
#
# The array task number is the study's index from `gwas2m-lncrna status`
# (run Step 1 / discovery first). Normally submitted for you by
# uq_submit_per_gwas.sh; by hand:
#
#     gwas2m-lncrna discover -p "coronary artery disease" -a EUR --max-studies 5
#     gwas2m-lncrna status   -p "coronary artery disease" -a EUR          # lists INDEX 1..N
#     PHENOTYPE="coronary artery disease" sbatch --array=1-5 examples/uq/uq_one_gwas.sh
#
# Defaults below are the UQ settings (partition "general", 24 h). A task that
# needs longer stops itself before the limit and resubmits just that task.
#
#SBATCH --job-name=gwas2m_gwas
#SBATCH --partition=general
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --output=gwas2m_gwas_%A_%a.log
##SBATCH --account=your_account

PHENOTYPE="${PHENOTYPE:-}"
if [[ -z "$PHENOTYPE" || -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo 'Submit as an array with the phenotype set, e.g.' >&2
    echo '  PHENOTYPE="coronary artery disease" sbatch --array=1-5 uq_one_gwas.sh' >&2
    exit 1
fi
ANCESTRY="${ANCESTRY:-EUR}"
# SLURM runs a spooled copy of this script, so ask the scheduler for the real path.
SCRIPT_PATH="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | awk -F= '/ Command=/{print $2; exit}' || true)"
SCRIPT_PATH="${SCRIPT_PATH:-$0}"
UQ_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
source "$UQ_DIR/uq_common.sh"

echo "Phenotype    : $PHENOTYPE ($ANCESTRY), study index $SLURM_ARRAY_TASK_ID"
run_with_resubmit gwas2m-lncrna run -p "$PHENOTYPE" -a "$ANCESTRY" --from-step 2 \
    --index "$SLURM_ARRAY_TASK_ID" --threads "${SLURM_CPUS_PER_TASK:-8}"

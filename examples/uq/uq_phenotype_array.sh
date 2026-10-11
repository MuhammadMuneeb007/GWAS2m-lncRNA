#!/bin/bash
# UQ HPC: one array task per phenotype listed in a text file (one per line).
#
#  1. Prepare reference data once:   sbatch examples/uq/uq_setup.sh
#  2. Then, after it finishes:
#        sbatch --array=1-$(grep -c . examples/phenotypes.txt) examples/uq/uq_phenotype_array.sh
#     (PHENOTYPE_FILE=my_list.txt sbatch --array=1-N ...  for your own list)
#  3. When every task has finished, write one combined report:
#        gwas2m-lncrna analysis -p "coronary artery disease" -p "Parkinson's disease" -a EUR \
#            --config examples/uq/gwas2m-lncrna.uq.yaml
#
# Defaults below are the UQ settings (partition "general", 24 h). A task that
# needs longer stops itself before the limit and resubmits just that task.
#
#SBATCH --job-name=gwas2m_array
#SBATCH --partition=general
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --output=gwas2m_array_%A_%a.log
##SBATCH --account=your_account

ANCESTRY="${ANCESTRY:-EUR}"
# SLURM runs a spooled copy of this script, so ask the scheduler for the real path.
SCRIPT_PATH="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | awk -F= '/ Command=/{print $2; exit}' || true)"
SCRIPT_PATH="${SCRIPT_PATH:-$0}"
UQ_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
source "$UQ_DIR/uq_common.sh"
PHENOTYPE_FILE="${PHENOTYPE_FILE:-$UQ_DIR/../phenotypes.txt}"

PHENOTYPE="$(grep . "$PHENOTYPE_FILE" | sed -n "${SLURM_ARRAY_TASK_ID}p")"
if [[ -z "$PHENOTYPE" ]]; then
    echo "No phenotype on line ${SLURM_ARRAY_TASK_ID} of $PHENOTYPE_FILE" >&2
    exit 1
fi
echo "Phenotype    : $PHENOTYPE ($ANCESTRY)"
run_with_resubmit gwas2m-lncrna run -p "$PHENOTYPE" -a "$ANCESTRY" --threads "${SLURM_CPUS_PER_TASK:-8}"

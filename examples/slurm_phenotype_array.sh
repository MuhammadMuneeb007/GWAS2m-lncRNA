#!/bin/bash
# One SLURM array task per phenotype listed in phenotypes.txt (one per line).
#
# 1. Prepare shared reference data once, so array tasks do not all download it:
#        gwas2m-lncrna setup --ancestry EUR
# 2. Submit with one task per line:
#        sbatch --array=1-$(grep -c . phenotypes.txt) slurm_phenotype_array.sh
#
#SBATCH --job-name=gwas2m-lncrna_array
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=gwas2m-lncrna_%A_%a.log
##SBATCH --partition=your_partition
##SBATCH --account=your_account

set -euo pipefail

ANCESTRY="EUR"
PHENOTYPE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" phenotypes.txt)"
if [[ -z "$PHENOTYPE" ]]; then
    echo "No phenotype on line ${SLURM_ARRAY_TASK_ID} of phenotypes.txt" >&2
    exit 1
fi

eval "$(mamba shell hook --shell bash 2>/dev/null || conda shell.bash hook)"
mamba activate gwas2m-lncrna 2>/dev/null || conda activate gwas2m-lncrna
cd "$SLURM_SUBMIT_DIR"

gwas2m-lncrna run -p "$PHENOTYPE" -a "$ANCESTRY" --threads "${SLURM_CPUS_PER_TASK:-4}"

# After every task has finished, build one combined report on a login node or in a small job:
#   gwas2m-lncrna analysis -p "coronary artery disease" -p "Parkinson's disease" -a EUR

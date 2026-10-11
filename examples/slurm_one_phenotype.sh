#!/bin/bash
# Run one phenotype end to end (steps 1-7) and write the report.
# Edit the #SBATCH lines for your cluster, then:  sbatch slurm_one_phenotype.sh
#
#SBATCH --job-name=gwas2m-lncrna
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=gwas2m-lncrna_%j.log
##SBATCH --partition=your_partition
##SBATCH --account=your_account

set -euo pipefail

PHENOTYPE="coronary artery disease"
ANCESTRY="EUR"

# Activate the environment created with: mamba env create -f environment.yml
eval "$(mamba shell hook --shell bash 2>/dev/null || conda shell.bash hook)"
mamba activate gwas2m-lncrna 2>/dev/null || conda activate gwas2m-lncrna

# Run from the folder holding gwas2m-lncrna.yaml (or add --config /path/to/gwas2m-lncrna.yaml).
cd "$SLURM_SUBMIT_DIR"

gwas2m-lncrna run -p "$PHENOTYPE" -a "$ANCESTRY" --threads "${SLURM_CPUS_PER_TASK:-4}" --analysis

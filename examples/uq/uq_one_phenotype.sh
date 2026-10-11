#!/bin/bash
# UQ HPC: run one phenotype end to end (steps 1-7) and write the report.
#
#     PHENOTYPE="coronary artery disease" sbatch examples/uq/uq_one_phenotype.sh
#     PHENOTYPE="Parkinson's disease" ANCESTRY=EUR sbatch examples/uq/uq_one_phenotype.sh
#     PHENOTYPE="coronary artery disease" MAX_STUDIES=5 sbatch examples/uq/uq_one_phenotype.sh   # top 5 GWAS
#
# Submit from the folder where results should go (or that holds gwas2m-lncrna.yaml).
# Defaults below are the UQ settings (partition "general", 24 h). If a run needs
# longer, the job stops itself before the limit and resubmits; finished studies
# and loci are skipped, so it continues where it stopped.
#
#SBATCH --job-name=gwas2m_lncrna
#SBATCH --partition=general
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --output=gwas2m_lncrna_%j.log
##SBATCH --account=your_account

PHENOTYPE="${PHENOTYPE:-}"
if [[ -z "$PHENOTYPE" ]]; then
    echo 'Set the phenotype, e.g.  PHENOTYPE="coronary artery disease" sbatch uq_one_phenotype.sh' >&2
    exit 1
fi
ANCESTRY="${ANCESTRY:-EUR}"
# SLURM runs a spooled copy of this script, so ask the scheduler for the real path.
SCRIPT_PATH="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | awk -F= '/ Command=/{print $2; exit}' || true)"
SCRIPT_PATH="${SCRIPT_PATH:-$0}"
UQ_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
source "$UQ_DIR/uq_common.sh"

echo "Phenotype    : $PHENOTYPE ($ANCESTRY)"
run_with_resubmit gwas2m-lncrna run -p "$PHENOTYPE" -a "$ANCESTRY" \
    --threads "${SLURM_CPUS_PER_TASK:-8}" --analysis ${DISCOVERY_ARGS[@]+"${DISCOVERY_ARGS[@]}"}

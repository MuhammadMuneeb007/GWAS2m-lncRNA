#!/bin/bash
# UQ HPC: download and build the shared reference data once
# (1000 Genomes LD panel, VEP 116 cache, GRCh38 genome, GENCODE).
#
#     sbatch examples/uq/uq_setup.sh            # European (EUR) panel
#     ANCESTRY=AFR sbatch examples/uq/uq_setup.sh
#
# Defaults below are the UQ settings (partition "general", 24 h). Override any
# of them on the command line, e.g.  sbatch --partition=other --time=12:00:00 ...
#
#SBATCH --job-name=gwas2m_setup
#SBATCH --partition=general
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=gwas2m_setup_%j.log
##SBATCH --account=your_account

ANCESTRY="${ANCESTRY:-EUR}"
# SLURM runs a spooled copy of this script, so ask the scheduler for the real path.
SCRIPT_PATH="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | awk -F= '/ Command=/{print $2; exit}' || true)"
SCRIPT_PATH="${SCRIPT_PATH:-$0}"
UQ_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
source "$UQ_DIR/uq_common.sh"

run_with_resubmit gwas2m-lncrna setup --ancestry "$ANCESTRY" --threads "${SLURM_CPUS_PER_TASK:-4}"

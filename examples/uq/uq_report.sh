#!/bin/bash
# UQ HPC: write the combined report (gwas2m-lncrna analysis) for the phenotypes
# listed in PHENOTYPE_FILE (one per line). Submitted by uq_submit_per_gwas.sh to
# run after the per-GWAS jobs; if per-GWAS jobs are still queued or running
# (e.g. tasks that resubmitted themselves), it waits for them and resubmits itself.
#
#     PHENOTYPE_FILE=examples/phenotypes.txt sbatch examples/uq/uq_report.sh
#
#SBATCH --job-name=gwas2m_report
#SBATCH --partition=general
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=gwas2m_report_%j.log
##SBATCH --account=your_account

ANCESTRY="${ANCESTRY:-EUR}"
SCRIPT_PATH="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | awk -F= '/ Command=/{print $2; exit}' || true)"
SCRIPT_PATH="${SCRIPT_PATH:-$0}"
UQ_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
source "$UQ_DIR/uq_common.sh"

PHENOTYPE_FILE="${PHENOTYPE_FILE:-$UQ_DIR/../phenotypes.txt}"

# Per-GWAS tasks may have resubmitted themselves after the 24 h limit; wait for them.
if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v squeue >/dev/null 2>&1; then
    WAITING="$( (squeue -h -u "$USER" -n gwas2m_gwas -o %A || true) | sort -u | paste -sd: -)"
    if [[ -n "$WAITING" ]]; then
        echo "Per-GWAS jobs still active ($WAITING); resubmitting the report to run after them."
        sbatch --dependency="afterany:$WAITING" --export=ALL "$SCRIPT_PATH"
        exit 0
    fi
fi

PHENOTYPE_ARGS=()
while IFS= read -r line; do
    [[ -n "${line// }" ]] && PHENOTYPE_ARGS+=(-p "$line")
done < "$PHENOTYPE_FILE"
echo "Report for   : $(( ${#PHENOTYPE_ARGS[@]} / 2 )) phenotypes from $PHENOTYPE_FILE ($ANCESTRY)"
gwas2m-lncrna analysis "${PHENOTYPE_ARGS[@]}" -a "$ANCESTRY"

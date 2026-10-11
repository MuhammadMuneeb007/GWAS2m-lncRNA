#!/bin/bash
# UQ HPC: run every GWAS in parallel - one SLURM array task per GWAS study.
# Run this on the LOGIN node (bash, not sbatch), from the folder where results go:
#
#     bash examples/uq/uq_submit_per_gwas.sh                         # phenotypes in examples/phenotypes.txt
#     MAX_STUDIES=5 bash examples/uq/uq_submit_per_gwas.sh           # top 5 GWAS per phenotype
#     PHENOTYPE_FILE=my_list.txt ANCESTRY=EUR bash examples/uq/uq_submit_per_gwas.sh
#
# What it does:
#   1. If the reference data is not ready, submits uq_setup.sh and makes every
#      other job wait for it.
#   2. For each phenotype: runs Step 1 (GWAS Catalog discovery) here, prints the
#      studies with their INDEX, and submits uq_one_gwas.sh as an array
#      --array=1-N (task i = study INDEX i, steps 2-7).
#   3. Submits uq_report.sh to write one combined report after all tasks finish.
#
# Options (environment variables):
#   MAX_STUDIES=N         top N GWAS per phenotype (else discovery.max_studies in the config)
#   REFRESH_DISCOVERY=1   search the GWAS Catalog again for phenotypes already discovered
#   MAX_PARALLEL=N        at most N GWAS jobs running at once per phenotype (array %N)
#   ANCESTRY=EUR          EUR, AFR, EAS, SAS or AMR
#   GWAS2M_LNCRNA_CONFIG  your configuration (else ./gwas2m-lncrna.yaml, else the UQ sample)

ANCESTRY="${ANCESTRY:-EUR}"
UQ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_PATH="${BASH_SOURCE[0]}"
source "$UQ_DIR/uq_common.sh"

PHENOTYPE_FILE="${PHENOTYPE_FILE:-$UQ_DIR/../phenotypes.txt}"
mapfile -t PHENOTYPES < <(grep -v '^[[:space:]]*$' "$PHENOTYPE_FILE")
if (( ${#PHENOTYPES[@]} == 0 )); then
    echo "No phenotypes in $PHENOTYPE_FILE" >&2
    exit 1
fi

# ---------------------------------------------------------------- 1. setup
DEPENDENCY=()
if ! gwas2m-lncrna check --ready --ancestry "$ANCESTRY"; then
    SETUP_JOB="$(ANCESTRY="$ANCESTRY" sbatch --parsable "$UQ_DIR/uq_setup.sh")"
    echo "Reference data missing: submitted setup job $SETUP_JOB; GWAS jobs will wait for it."
    DEPENDENCY=(--dependency="afterok:$SETUP_JOB")
fi

DISCOVER_ARGS=()
[[ -n "${MAX_STUDIES:-}" ]] && DISCOVER_ARGS+=(--max-studies "$MAX_STUDIES")
[[ "${REFRESH_DISCOVERY:-0}" == 1 ]] && DISCOVER_ARGS+=(--force)
THROTTLE=""
[[ -n "${MAX_PARALLEL:-}" ]] && THROTTLE="%$MAX_PARALLEL"

# ------------------------------------------------- 2. discovery + arrays
ARRAY_JOBS=()
SUBMITTED=()
for PHENOTYPE in "${PHENOTYPES[@]}"; do
    echo
    echo "=================================================================="
    echo "$PHENOTYPE ($ANCESTRY)"
    echo "=================================================================="
    if ! gwas2m-lncrna discover -p "$PHENOTYPE" -a "$ANCESTRY" ${DISCOVER_ARGS[@]+"${DISCOVER_ARGS[@]}"}; then
        echo "Discovery failed for '$PHENOTYPE'; skipping it." >&2
        continue
    fi
    N="$(gwas2m-lncrna status -p "$PHENOTYPE" -a "$ANCESTRY" --count)"
    gwas2m-lncrna status -p "$PHENOTYPE" -a "$ANCESTRY"
    JOB="$(PHENOTYPE="$PHENOTYPE" ANCESTRY="$ANCESTRY" sbatch --parsable ${DEPENDENCY[@]+"${DEPENDENCY[@]}"} \
           --array="1-${N}${THROTTLE}" "$UQ_DIR/uq_one_gwas.sh")"
    echo "Submitted $N GWAS jobs for '$PHENOTYPE': array job $JOB (task i = INDEX i above)"
    ARRAY_JOBS+=("${JOB%%;*}")
    SUBMITTED+=("$PHENOTYPE")
done

if (( ${#ARRAY_JOBS[@]} == 0 )); then
    echo "Nothing was submitted." >&2
    exit 1
fi

# --------------------------------------------------------------- 3. report
LIST="$PWD/gwas2m_report_phenotypes_$(date +%Y%m%d_%H%M%S).txt"
printf '%s\n' "${SUBMITTED[@]}" > "$LIST"
REPORT_JOB="$(PHENOTYPE_FILE="$LIST" ANCESTRY="$ANCESTRY" sbatch --parsable \
              --dependency="afterany:$(IFS=:; echo "${ARRAY_JOBS[*]}")" "$UQ_DIR/uq_report.sh")"
echo
echo "Report job $REPORT_JOB runs after all GWAS jobs (phenotypes: $LIST)."
echo "Follow progress with:  squeue -u \$USER   and   gwas2m-lncrna status -p \"<phenotype>\" -a $ANCESTRY"

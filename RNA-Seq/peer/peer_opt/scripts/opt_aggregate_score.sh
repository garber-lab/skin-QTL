#!/bin/bash
# =============================================================================
# opt_aggregate_score.sh  --  after all sub-chunks of one (CT, setting, fold)
# have been fit, run QC + eigenMT per sub-chunk, aggregate to fold level, then
# BH + score. QC/eigenMT are cheap vs fitting, so they run here in one job.
#
# Args: CT METHOD K nGPC FOLD COV_DIR FOLD_DIR GENELIST_DIR RUN_DIR SUMMARY SCRIPTDIR SIF
# Env: CONDITIONS (default "PBS"), QTLTYPES (default "eQTL")
# =============================================================================
set -uo pipefail
CT=$1; METHOD=$2; K=$3; nGPC=$4; FOLD=$5
COV_DIR=$6; FOLD_DIR=$7; GENELIST_DIR=$8; RUN_DIR=$9; SUMMARY=${10}; SCRIPTDIR=${11}; SIF=${12}
CONDITIONS=${CONDITIONS:-"PBS"}; QTLTYPES=${QTLTYPES:-"eQTL"}

Kpad=$(printf "%02d" "$K")
SETTING_ID="${METHOD}_K${Kpad}_gpc${nGPC}"
BASE="${RUN_DIR}/${SETTING_ID}/${FOLD}"
mkdir -p "${BASE}/${CT}/results_QC" "${BASE}/${CT}/eigenMT/results" "${BASE}/results"

# Singularity/Apptainer only forwards host vars carrying the SINGULARITYENV_ /
# APPTAINERENV_ prefix into the container. Set the plain names too (harmless).
export OPT_ROOTDIR="$BASE" OPT_CHUNKS_DIR="$FOLD_DIR" OPT_GENELIST_DIR="$GENELIST_DIR"
export SINGULARITYENV_OPT_ROOTDIR="$BASE" SINGULARITYENV_OPT_CHUNKS_DIR="$FOLD_DIR" SINGULARITYENV_OPT_GENELIST_DIR="$GENELIST_DIR"
export APPTAINERENV_OPT_ROOTDIR="$BASE" APPTAINERENV_OPT_CHUNKS_DIR="$FOLD_DIR" APPTAINERENV_OPT_GENELIST_DIR="$GENELIST_DIR"
R() { singularity exec "$SIF" Rscript "$@"; }

mapfile -t SUBS < "${FOLD_DIR}/subchunks_${FOLD}.txt"
N_GENES=$(wc -l < "${FOLD_DIR}/fold_genes_${FOLD}.txt")
STATUS="ok"; t0=$(date +%s)

for cond in $CONDITIONS; do
  for qt in $QTLTYPES; do
    [[ "$qt" == "reQTL" && "$cond" == "PBS" ]] && continue
    mkdir -p "${BASE}/${CT}/eigenMT/${cond}/${qt}"   # step3 writes here but doesn't create it

    # per sub-chunk: QC then eigenMT (chunk_id = sub tag)
    for sub in "${SUBS[@]}"; do
      [[ -f "${BASE}/${CT}/results/result_${sub}.tsv" ]] || { STATUS="missing_fit_${sub}"; continue; }
      R "${SCRIPTDIR}/step2_QTL_result_QC.R" "$CT" "$cond" "$qt" "$sub" || STATUS="qc_failed_${sub}"
      R "${SCRIPTDIR}/step3_QTL_result_correction_eigenMT.R" "$CT" "$cond" "$qt" "$sub" || STATUS="mt_failed_${sub}"
    done

    # aggregate sub-chunk QC -> fold-level QC (for lambda), header kept once
    QCOUT="${BASE}/${CT}/results_QC/modeling_stats_postQC_${CT}_${cond}_${qt}_${FOLD}.txt"
    QCPARTS=( "${BASE}/${CT}/results_QC/modeling_stats_postQC_${CT}_${cond}_${qt}_${FOLD}"_s*.txt )
    awk 'FNR==1{if(seen++)next}1' "${QCPARTS[@]}" > "$QCOUT" 2>/dev/null || STATUS="qc_aggregate_failed"

    # aggregate sub-chunk eigenMT -> fold-level eigenMT (step4's expected name)
    MTOUT="${BASE}/${CT}/eigenMT/results/${CT}_${cond}_${qt}.eigenMT.txt"
    MTPARTS=( "${BASE}/${CT}/eigenMT/${cond}/${qt}/${CT}_${cond}_${qt}_${FOLD}"_s*.eigenMT.tsv )
    if awk 'FNR==1{if(seen++)next}1' "${MTPARTS[@]}" > "$MTOUT" 2>/dev/null; then
      R "${SCRIPTDIR}/step4_QTL_result_correction_BH.R" "${CT}_${cond}_${qt}.eigenMT.txt" || STATUS="bh_failed_${cond}_${qt}"
    else
      STATUS="no_eigenmt_${cond}_${qt}"
    fi
  done
done

t1=$(date +%s); RUNTIME=$(( t1 - t0 ))
SEED_FILE="${FOLD_DIR}/fold_seed_egenes_${FOLD}.txt"
SEED_ARG=(); [[ -f "$SEED_FILE" ]] && SEED_ARG=(--seed-genes-file "$SEED_FILE")

echo "[$(date)] SCORE $CT $SETTING_ID $FOLD (status=$STATUS, agg ${RUNTIME}s)"
R "${SCRIPTDIR}/opt_lambda_and_summary.R" --base-dir "$BASE" --cell-type "$CT" --fold "$FOLD" \
  --method "$METHOD" --K "$K" --nGPC "$nGPC" --summary-file "$SUMMARY" \
  --runtime-sec "$RUNTIME" --status "$STATUS" --n-genes-tested "$N_GENES" \
  "${SEED_ARG[@]}"
echo "[$(date)] DONE $CT $SETTING_ID $FOLD"

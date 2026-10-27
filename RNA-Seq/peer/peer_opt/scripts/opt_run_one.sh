#!/bin/bash
# =============================================================================
# opt_run_one.sh  --  run the full covariate-optimization pipeline for ONE
#                     (cell type, covariate setting, fold).
#
# Chains: opt_fit_qtl_subset.R -> step2 (QC) -> step3 (eigenMT) -> step4 (BH)
#         -> opt_lambda_and_summary.R, all under a per-(setting,fold) base-dir
#         so nothing collides with production or other settings.
#
# Reuses your real steps 2-4 via the OPT_ROOTDIR / OPT_CHUNKS_DIR /
# OPT_GENELIST_DIR env hooks (unset = production; here we set them).
#
# Args (positional):
#   1 CT           FRB|KRT|MEL
#   2 METHOD       peer | expression_pc | expression_pc_condremoved
#   3 K            e.g. 8
#   4 nGPC         e.g. 2
#   5 FOLD         e.g. fold1
#   6 COV_DIR      covariate output dir (…/peer_opt) from generate_covariates.R
#   7 FOLD_DIR     dir with per-fold pairs/genotype (from opt_make_folds.R)
#   8 GENELIST_DIR dir with gene_list_<CT>_<cond>_{expressed,induced}.txt
#   9 RUN_DIR      root for optimizer outputs
#  10 SUMMARY      master summary TSV to append to
#  11 SCRIPTDIR    dir containing opt_*.R and step2/3/4 R scripts
#  12 SIF          singularity image
# Optional env: CONDITIONS (default "PBS IFNG IFNB TNF"), QTLTYPES (default "eQTL reQTL"),
#               NWORKERS (default 5)
# =============================================================================
set -euo pipefail

CT=$1; METHOD=$2; K=$3; nGPC=$4; FOLD=$5
COV_DIR=$6; FOLD_DIR=$7; GENELIST_DIR=$8; RUN_DIR=$9; SUMMARY=${10}
SCRIPTDIR=${11}; SIF=${12}
CONDITIONS=${CONDITIONS:-"PBS IFNG IFNB TNF"}
QTLTYPES=${QTLTYPES:-"eQTL reQTL"}
NWORKERS=${NWORKERS:-5}

Kpad=$(printf "%02d" "$K")
case "$METHOD" in
  peer)                       SUBDIR=peer;                         PREFIX=peer   ;;
  expression_pc)              SUBDIR=expression_pcs;               PREFIX=exprPC ;;
  expression_pc_condremoved)  SUBDIR=expression_pcs_condremoved;   PREFIX=exprPC ;;
  *) echo "Unknown method: $METHOD" >&2; exit 1 ;;
esac
COV_FILE="${COV_DIR}/covariates/${SUBDIR}/${CT}/${PREFIX}_K${Kpad}.tsv"
[[ -f "$COV_FILE" ]] || { echo "Covariate file missing: $COV_FILE" >&2; exit 1; }

SETTING_ID="${METHOD}_K${Kpad}_gpc${nGPC}"
BASE="${RUN_DIR}/${SETTING_ID}/${FOLD}"      # per-(setting,fold) pipeline root
mkdir -p "${BASE}/${CT}/results" "${BASE}/${CT}/results_QC" \
         "${BASE}/${CT}/eigenMT/results" "${BASE}/results"

PAIRS="${FOLD_DIR}/pairs_chunk_${FOLD}.tsv"
GENO="${FOLD_DIR}/genotype_pairs_chunk_${FOLD}.tsv"
RESULT="${BASE}/${CT}/results/result_${FOLD}.tsv"
N_GENES=$(( $(wc -l < "${FOLD_DIR}/fold_genes_${FOLD}.txt") ))
N_PAIRS=$(( $(wc -l < "$PAIRS") - 1 ))

# env hooks so steps 2-4 write under BASE and read per-fold chunks + shared gene lists
export OPT_ROOTDIR="$BASE"
export OPT_CHUNKS_DIR="$FOLD_DIR"
export OPT_GENELIST_DIR="$GENELIST_DIR"

R() { singularity exec "$SIF" Rscript "$@"; }
STATUS="ok"; t0=$(date +%s)

echo "[$(date)] FIT  $CT $SETTING_ID $FOLD"
R "${SCRIPTDIR}/opt_fit_qtl_subset.R" "$CT" "$COV_FILE" "$nGPC" "$PAIRS" "$GENO" "$RESULT" "$NWORKERS" \
  || STATUS="fit_failed"

if [[ "$STATUS" == "ok" ]]; then
  for cond in $CONDITIONS; do
    for qt in $QTLTYPES; do
      # reQTL is only defined for stimulated conditions
      if [[ "$qt" == "reQTL" && "$cond" == "PBS" ]]; then continue; fi
      echo "[$(date)] QC   $CT $cond $qt $FOLD"
      R "${SCRIPTDIR}/step2_QTL_result_QC.R" "$CT" "$cond" "$qt" "$FOLD" || { STATUS="qc_failed_${cond}_${qt}"; continue; }

      echo "[$(date)] MT   $CT $cond $qt $FOLD"
      R "${SCRIPTDIR}/step3_QTL_result_correction_eigenMT.R" "$CT" "$cond" "$qt" "$FOLD" \
        || { STATUS="eigenmt_failed_${cond}_${qt}"; continue; }

      # bridge step3 per-fold output -> step4's expected aggregated eigenMT file
      SRC="${BASE}/${CT}/eigenMT/${cond}/${qt}/${CT}_${cond}_${qt}_${FOLD}.eigenMT.tsv"
      DST="${BASE}/${CT}/eigenMT/results/${CT}_${cond}_${qt}.eigenMT.txt"
      if [[ -f "$SRC" ]]; then cp "$SRC" "$DST"; else STATUS="no_eigenmt_${cond}_${qt}"; continue; fi

      echo "[$(date)] BH   $CT $cond $qt $FOLD"
      R "${SCRIPTDIR}/step4_QTL_result_correction_BH.R" "${CT}_${cond}_${qt}.eigenMT.txt" \
        || STATUS="bh_failed_${cond}_${qt}"
    done
  done
fi

t1=$(date +%s); RUNTIME=$(( t1 - t0 ))

echo "[$(date)] SCORE $CT $SETTING_ID $FOLD (status=$STATUS, ${RUNTIME}s)"
R "${SCRIPTDIR}/opt_lambda_and_summary.R" --base-dir "$BASE" --cell-type "$CT" --fold "$FOLD" \
  --method "$METHOD" --K "$K" --nGPC "$nGPC" --summary-file "$SUMMARY" \
  --runtime-sec "$RUNTIME" --status "$STATUS" --n-genes-tested "$N_GENES" --n-pairs "$N_PAIRS"

echo "[$(date)] DONE $CT $SETTING_ID $FOLD"

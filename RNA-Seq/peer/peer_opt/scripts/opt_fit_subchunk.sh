#!/bin/bash
# =============================================================================
# opt_fit_subchunk.sh  --  fit ONE sub-chunk of one (CT, setting, fold).
# Each sub-chunk is production-sized (~50 genes), so geno tables are small and
# this runs like a normal chunk job. Aggregation happens later per fold.
#
# Args: CT METHOD K nGPC FOLD SUBCHUNK COV_DIR FOLD_DIR RUN_DIR SCRIPTDIR SIF
# Optional env: NWORKERS (default 5)
# =============================================================================
set -euo pipefail
CT=$1; METHOD=$2; K=$3; nGPC=$4; FOLD=$5; SUB=$6
COV_DIR=$7; FOLD_DIR=$8; RUN_DIR=$9; SCRIPTDIR=${10}; SIF=${11}
NWORKERS=${NWORKERS:-5}

Kpad=$(printf "%02d" "$K")
case "$METHOD" in
  peer)                      SUBDIR=peer;                       PREFIX=peer   ;;
  expression_pc)             SUBDIR=expression_pcs;             PREFIX=exprPC ;;
  expression_pc_condremoved) SUBDIR=expression_pcs_condremoved; PREFIX=exprPC ;;
  *) echo "Unknown method: $METHOD" >&2; exit 1 ;;
esac
COV_FILE="${COV_DIR}/covariates/${SUBDIR}/${CT}/${PREFIX}_K${Kpad}.tsv"
[[ -f "$COV_FILE" ]] || { echo "Covariate file missing: $COV_FILE" >&2; exit 1; }

SETTING_ID="${METHOD}_K${Kpad}_gpc${nGPC}"
BASE="${RUN_DIR}/${SETTING_ID}/${FOLD}"
mkdir -p "${BASE}/${CT}/results"

PAIRS="${FOLD_DIR}/pairs_chunk_${SUB}.tsv"
GENO="${FOLD_DIR}/genotype_pairs_chunk_${SUB}.tsv"
RESULT="${BASE}/${CT}/results/result_${SUB}.tsv"

if [[ -f "$RESULT" ]]; then echo "exists, skip: $RESULT"; exit 0; fi   # idempotent rerun
echo "[$(date)] FIT $CT $SETTING_ID $SUB"
singularity exec "$SIF" Rscript "${SCRIPTDIR}/opt_fit_qtl_subset.R" \
  "$CT" "$COV_FILE" "$nGPC" "$PAIRS" "$GENO" "$RESULT" "$NWORKERS"

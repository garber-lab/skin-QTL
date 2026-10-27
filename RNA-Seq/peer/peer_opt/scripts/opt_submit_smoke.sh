#!/bin/bash
# =============================================================================
# opt_submit.sh  --  two-phase covariate-optimization sweep.
#   Phase FIT: array over (CT, method, K, nGPC, fold, sub-chunk) -- the heavy,
#              parallel part (each sub-chunk ~ one production chunk).
#   Phase AGG: array over (CT, method, K, nGPC, fold), dependent on FIT -- runs
#              QC + eigenMT per sub-chunk, aggregates to fold level, BH + score.
#
# STEP 0 (once): build PER-CELL-TYPE seeded, sub-chunked folds (see below).
#   PROD=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models
#   for ct in FRB KRT MEL; do
#     singularity exec $SIF Rscript opt_make_folds.R \
#       --chunks-dir $PROD/data/chunks --out-dir ${DIR}/peer_opt/folds/$ct \
#       --n-folds 3 --n-genes-per-fold 1000 \
#       --seed-genes-file $PROD/results/${ct}_PBS_eQTL_gene_fdr05_genelist.txt \
#       --n-seed-per-fold 100 --subchunk-genes 50 --seed 42
#   done
#
# Then:  bash opt_submit.sh          (SMOKE=1 bash opt_submit.sh to test one fold)
# =============================================================================
set -euo pipefail

## ------------------------------ CONFIG ------------------------------------ ##
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer
SCRIPTDIR=${DIR}/peer_opt/scripts
SIF=/share/pkg/containers/rstudio_example/r_4.5.2.sif

COV_DIR=${DIR}/peer_opt
FOLD_BASE=${DIR}/peer_opt/folds_smoke
GENELIST_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/data
RUN_DIR=${DIR}/peer_opt/optimizer
SUMMARY=${RUN_DIR}/covariate_optimization_summary.tsv

export CONDITIONS="PBS"          # Option 2: PBS eQTL only
export QTLTYPES="eQTL"

CTS="FRB"
METHODS="expression_pc_condremoved"
KGRID="5"
NGPC="3"
FOLDS="fold1"

FIT_PARALLEL=120                  # like your production array
SMOKE=${SMOKE:-0}
## -------------------------------------------------------------------------- ##

mkdir -p "$RUN_DIR"
FITC="${RUN_DIR}/fit_combos.tsv"; AGGC="${RUN_DIR}/agg_combos.tsv"

build_combos() {
  : > "$FITC"; : > "$AGGC"
  for ct in $CTS; do for m in $METHODS; do for k in $KGRID; do for g in $NGPC; do for f in $FOLDS; do
    printf "%s\t%s\t%s\t%s\t%s\n" "$ct" "$m" "$k" "$g" "$f" >> "$AGGC"
    subfile="${FOLD_BASE}/${ct}/subchunks_${f}.txt"
    [[ -f "$subfile" ]] || { echo "Missing $subfile -- run opt_make_folds.R first" >&2; exit 1; }
    while read -r sub; do
      [[ -n "$sub" ]] && printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$ct" "$m" "$k" "$g" "$f" "$sub" >> "$FITC"
    done < "$subfile"
  done; done; done; done; done
}

run_fit_line() { local l; l=$(sed -n "${1}p" "$FITC"); read -r ct m k g f sub <<< "$l"
  bash "${SCRIPTDIR}/opt_fit_subchunk.sh" "$ct" "$m" "$k" "$g" "$f" "$sub" \
    "$COV_DIR" "${FOLD_BASE}/${ct}" "$RUN_DIR" "$SCRIPTDIR" "$SIF"; }

run_agg_line() { local l; l=$(sed -n "${1}p" "$AGGC"); read -r ct m k g f <<< "$l"
  bash "${SCRIPTDIR}/opt_aggregate_score.sh" "$ct" "$m" "$k" "$g" "$f" \
    "$COV_DIR" "${FOLD_BASE}/${ct}" "$GENELIST_DIR" "$RUN_DIR" "$SUMMARY" "$SCRIPTDIR" "$SIF"; }

# Array-element mode (dispatched by bsub): run one line of the phase's combo file.
if [[ -n "${LSB_JOBINDEX:-}" ]]; then
  sleep $((RANDOM % 20))
  case "${PHASE:-}" in
    fit) run_fit_line "$LSB_JOBINDEX" ;;
    agg) run_agg_line "$LSB_JOBINDEX" ;;
    *) echo "PHASE unset" >&2; exit 1 ;;
  esac
  exit 0
fi

build_combos
NFIT=$(wc -l < "$FITC"); NAGG=$(wc -l < "$AGGC")
echo "Built $NFIT fit sub-chunk jobs, $NAGG aggregate jobs"

if [[ "$SMOKE" == "1" ]]; then
  echo "SMOKE: fitting all sub-chunks of the first fold, then aggregating ..."
  read -r ct m k g f <<< "$(sed -n '1p' "$AGGC")"
  while read -r sub; do
    [[ -n "$sub" ]] && bash "${SCRIPTDIR}/opt_fit_subchunk.sh" "$ct" "$m" "$k" "$g" "$f" "$sub" \
      "$COV_DIR" "${FOLD_BASE}/${ct}" "$RUN_DIR" "$SCRIPTDIR" "$SIF"
  done < "${FOLD_BASE}/${ct}/subchunks_${f}.txt"
  bash "${SCRIPTDIR}/opt_aggregate_score.sh" "$ct" "$m" "$k" "$g" "$f" \
    "$COV_DIR" "${FOLD_BASE}/${ct}" "$GENELIST_DIR" "$RUN_DIR" "$SUMMARY" "$SCRIPTDIR" "$SIF"
  echo "SMOKE done. Inspect: $SUMMARY"
  exit 0
fi

LOGDIR=${RUN_DIR}/logs; mkdir -p "$LOGDIR"
bsub -J "covfit[1-${NFIT}]%${FIT_PARALLEL}" -n "${NWORKERS:-5}" -R "span[hosts=1]" \
     -R "rusage[mem=4000]" -W 12:00 -q long \
     -o "${LOGDIR}/covfit_%J_%I.out" -e "${LOGDIR}/covfit_%J_%I.err" \
     "PHASE=fit bash ${BASH_SOURCE[0]}"
bsub -J "covagg[1-${NAGG}]%50" -w "done(covfit)" -n 2 -R "span[hosts=1]" \
     -R "rusage[mem=4000]" -W 8:00 -q long \
     -o "${LOGDIR}/covagg_%J_%I.out" -e "${LOGDIR}/covagg_%J_%I.err" \
     "PHASE=agg bash ${BASH_SOURCE[0]}"
echo "Submitted covfit[1-${NFIT}]%${FIT_PARALLEL} and covagg[1-${NAGG}] (waits on covfit)."
echo "Summary accrues in: $SUMMARY"

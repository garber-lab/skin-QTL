#!/bin/bash
# Contribution scores across ALL peaks, for every fold of every model.
# Prerequisite for any genome-wide motif scan: Fi-NeMo calls hits on contributions, and the
# pipeline only computed them for a 30K-peak subsample (chrombpnet_model/auxiliary/interpret_subsample).
#
# Usage:  bash step7_motifScan_contribs.sh [models...]        (default: all six)
#         DRYRUN=1 bash step7_motifScan_contribs.sh           (write job scripts, submit nothing)
#
# Output: ${CB}/motifscan/contribs/<model>/<model>_fold<f>.{counts,profile}_scores.h5 (+ bigwigs)
set -euo pipefail

DRYRUN=${DRYRUN:-0}
CB=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet

############################ CONFIG ############################
MODELS_DEFAULT="FRB_IFN FRB_PBS KRT_IFN KRT_PBS MEL_IFN MEL_PBS"
FOLDS="0 1 2 3 4"
BIAS=0.5
HEADS="counts profile"                               # -pc for chrombpnet contribs_bw ("counts" alone
                                                     # halves the cost; counts drives accessibility)
# Peaks to score. Default: that model's full peak set (open regions only, never the whole genome).
# Set PEAKS_OVERRIDE to a 10-column narrowPeak subset (e.g. peaks near your variants) to cut cost:
# runtime scales linearly with the number of regions.
PEAKS_OVERRIDE=${PEAKS_OVERRIDE:-}
OUT=${CB}/motifscan/contribs
DL=${CB}/ATAC_KRT_IFN_merged_alldonors/biasFactor${BIAS}/data/downloads
GENOME=${DL}/hg38.fa
CHROMSIZES=${DL}/hg38.chrom.sizes
SIF=${CB}/chrombpnet_latest.sif
CONDA_SITE=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages
HOST_SETUP="module load conda; source activate fastQTL"
RUN="singularity exec --nv --env PYTHONPATH=${CONDA_SITE} --env TF_FORCE_GPU_ALLOW_GROWTH=true ${SIF}"
# DeepSHAP over every peak (~150-250K regions) is the long pole: hours per fold on a V100
BSUB_OPTS='-q gpu -W 48:00 -n 4 -R "rusage[mem=16000]" -R "span[hosts=1]"'
################################################################

MODELS=${@:-${MODELS_DEFAULT}}
mkdir -p "${OUT}" "${CB}/motifscan/logs" "${CB}/motifscan/jobs"

for m in ${MODELS}; do
  base=${CB}/ATAC_${m}_merged_alldonors/biasFactor${BIAS}
  peaks=${PEAKS_OVERRIDE:-${base}/data/downloads/peaks_no_blacklist.narrowPeak}   # 10-column narrowPeak
  [[ -s "${peaks}" ]] || { echo "ERROR: peaks not found: ${peaks}" >&2; exit 1; }
  echo "${m}: $(wc -l < ${peaks}) peaks"

  jobF=${CB}/motifscan/jobs/contribs_${m}.sh
  cat > "${jobF}" <<EOF
#!/bin/bash
${HOST_SETUP}
set -euo pipefail
shopt -s nullglob
mkdir -p ${OUT}/${m}
for f in ${FOLDS}; do
  hits=(${base}/fold\${f}/chrombpnet_model/models/*chrombpnet_nobias.h5)
  if [[ \${#hits[@]} -ne 1 ]]; then
    echo "ERROR: expected 1 *chrombpnet_nobias.h5 in fold\${f}, found \${#hits[@]}" >&2; exit 1
  fi
  op=${OUT}/${m}/${m}_fold\${f}
  if [[ -s \${op}.counts_scores.h5 && -s \${op}.profile_scores.h5 ]]; then
    echo "fold\${f}: already done, skipping"; continue
  fi
  echo "== fold\${f} \${hits[0]}"; date
  ${RUN} bash -c "${ENV_SETUP:-} chrombpnet contribs_bw \
      -m \${hits[0]} -r ${peaks} -g ${GENOME} -c ${CHROMSIZES} \
      -op \${op} -pc ${HEADS}"
  date
done
echo "done ${m}"
EOF
  if [[ ${DRYRUN} -eq 1 ]]; then echo "[dry run] wrote ${jobF}"; continue; fi
  eval bsub ${BSUB_OPTS} -J "mscan_contribs_${m}" \
    -oo "${CB}/motifscan/logs/contribs_${m}.out" -eo "${CB}/motifscan/logs/contribs_${m}.err" \
    bash "${jobF}"
done

echo "Contribution scores -> ${OUT}/<model>/<model>_fold<f>.{counts,profile}_scores.h5"
echo "Jobs resume where they left off: reruns skip folds whose h5 files already exist."

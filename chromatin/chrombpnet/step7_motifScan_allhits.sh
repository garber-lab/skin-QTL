#!/bin/bash
# Call hits for ALL MoDISco motifs across the scored peaks, per model (one Fi-NeMo run each).
# Calling every motif costs about the same as calling one, so this replaces the per-TF scan:
# filter the resulting hits table by motif_name for any TF or pair you care about.
#
# Requires contribution scores from step7_motifScan_contribs.sh.
#
# Usage:  bash step7_motifScan_allhits.sh [models...]        (default: all six)
#         DRYRUN=1 bash step7_motifScan_allhits.sh
#
# Output per model, under ${CB}/motifscan/allhits_<head>/<model>/:
#   regions.npz            peak regions + fold-averaged contributions (cached, reused by later runs)
#   hits.tsv               every motif instance, with genomic coordinates and readable names
#   hits.bed               chr, start, end, motif, coefficient, strand (sorted)
#   motif_counts.tsv       instances per motif
set -euo pipefail

DRYRUN=${DRYRUN:-0}
CB=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet

############################ CONFIG ############################
MODELS_DEFAULT="FRB_IFN FRB_PBS KRT_IFN KRT_PBS MEL_IFN MEL_PBS"
BIAS=0.5
HEAD=profile                                         # must match the MoDISco report and contributions
CONTRIBS=${CB}/motifscan/contribs                    # from step7_motifScan_contribs.sh
OUT=${CB}/motifscan/allhits_${HEAD}
NAMES_PY=${CB}/step7_motifScan_names.py
NAMES_DIR=${CB}/motifscan/names_${HEAD}
PEAKS_OVERRIDE=${PEAKS_OVERRIDE:-}                   # same subset used for the contributions, if any
REGION_WIDTH=1000
FINEMO=/home/shuo.shan-umw/.conda/envs/finemo/bin/finemo
FINEMO_DEVICE=gpu                                    # many regions: GPU is worth it (torch cu126)
FINEMO_OPTS="-c 0.0001 -S 10000"
MODISCO_GLOB="fold0/chrombpnet_model/auxiliary/interpret_subsample/*modisco_results_${HEAD}_scores.h5"
SIF=${CB}/chrombpnet_latest.sif
CONDA_SITE=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages
HOST_SETUP="module load conda; source activate fastQTL"
RUN="singularity exec --env PYTHONPATH=${CONDA_SITE} ${SIF}"
BSUB_OPTS='-q gpu -W 8:00 -n 4 -R "rusage[mem=64000]" -R "span[hosts=1]"'
################################################################

MODELS=${@:-${MODELS_DEFAULT}}
[[ -f "${NAMES_PY}" ]] || { echo "ERROR: missing ${NAMES_PY}" >&2; exit 1; }
[[ -x "${FINEMO}" ]] || { echo "ERROR: finemo not executable: ${FINEMO}" >&2; exit 1; }
if [[ "${FINEMO_DEVICE}" == "cpu" ]]; then FINEMO_ENV='CUDA_VISIBLE_DEVICES=""'; else FINEMO_ENV=""; fi
mkdir -p "${OUT}" "${CB}/motifscan/logs" "${CB}/motifscan/jobs"

### readable motif names for every pattern of every model (seconds, runs here not in the job)
${RUN} bash -c "export PYTHONPATH=${CONDA_SITE}; python ${NAMES_PY} \
    --cb ${CB} --bias ${BIAS} --models ${MODELS// /,} --head ${HEAD} --out ${NAMES_DIR}"

for m in ${MODELS}; do
  base=${CB}/ATAC_${m}_merged_alldonors/biasFactor${BIAS}
  peaks=${PEAKS_OVERRIDE:-${base}/data/downloads/peaks_no_blacklist.narrowPeak}
  jobF=${CB}/motifscan/jobs/allhits_${HEAD}_${m}.sh
  cat > "${jobF}" <<EOF
#!/bin/bash
${HOST_SETUP}
set -euo pipefail
shopt -s nullglob
d=${OUT}/${m}
mkdir -p \${d}

mh=(${base}/${MODISCO_GLOB})
if [[ \${#mh[@]} -ne 1 ]]; then echo "ERROR: expected 1 MoDISco h5 for ${m}, found \${#mh[@]}" >&2; exit 1; fi
h5s=(${CONTRIBS}/${m}/${m}_fold*.${HEAD}_scores.h5)
if [[ \${#h5s[@]} -eq 0 ]]; then
  echo "ERROR: no ${HEAD} contributions for ${m} (run step7_motifScan_contribs.sh first)" >&2; exit 1
fi
echo "${m}: \${#h5s[@]} fold(s), \$(wc -l < ${peaks}) peaks"; date

# peak regions with genomic coordinates, contributions averaged across folds (cached)
if [[ ! -s \${d}/regions.npz ]]; then
  ${FINEMO_ENV} ${FINEMO} extract-regions-chrombpnet-h5 \
    -c "\${h5s[@]}" -p ${peaks} -w ${REGION_WIDTH} -o \${d}/regions.npz
fi
date

# every motif at once (no -I): same cost as calling a single motif
${FINEMO_ENV} ${FINEMO} call-hits -r \${d}/regions.npz -m \${mh[0]} -o \${d} \
  -N ${NAMES_DIR}/${m}/names.txt ${FINEMO_OPTS}
echo "\${mh[0]}" > \${d}/modisco_source.txt
date

# BED of instances + per-motif counts
awk -F'\t' 'NR>1 && \$1!="NA" {print \$1"\t"\$2"\t"\$3"\t"\$6"\t"\$7"\t"\$13}' \${d}/hits.tsv \
  | sort -k1,1 -k2,2n > \${d}/hits.bed
awk -F'\t' 'NR>1 {c[\$6]++} END {for (k in c) print k"\t"c[k]}' \${d}/hits.tsv \
  | sort -k2,2nr > \${d}/motif_counts.tsv
echo "${m}: \$(wc -l < \${d}/hits.bed) instances, \$(wc -l < \${d}/motif_counts.tsv) motifs"
head -15 \${d}/motif_counts.tsv
EOF
  if [[ ${DRYRUN} -eq 1 ]]; then echo "[dry run] wrote ${jobF}"; continue; fi
  eval bsub ${BSUB_OPTS} -J "mscan_allhits_${m}" \
    -oo "${CB}/motifscan/logs/allhits_${HEAD}_${m}.out" -eo "${CB}/motifscan/logs/allhits_${HEAD}_${m}.err" \
    bash "${jobF}"
done

echo "All-motif hits -> ${OUT}/<model>/hits.tsv (+ hits.bed, motif_counts.tsv)"

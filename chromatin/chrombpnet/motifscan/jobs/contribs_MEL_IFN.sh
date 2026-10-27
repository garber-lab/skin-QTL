#!/bin/bash
module load conda; source activate fastQTL
set -euo pipefail
shopt -s nullglob
mkdir -p /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/motifscan/contribs/MEL_IFN
for f in 0 1 2 3 4; do
  hits=(/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5/fold${f}/chrombpnet_model/models/*chrombpnet_nobias.h5)
  if [[ ${#hits[@]} -ne 1 ]]; then
    echo "ERROR: expected 1 *chrombpnet_nobias.h5 in fold${f}, found ${#hits[@]}" >&2; exit 1
  fi
  op=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/motifscan/contribs/MEL_IFN/MEL_IFN_fold${f}
  if [[ -s ${op}.counts_scores.h5 && -s ${op}.profile_scores.h5 ]]; then
    echo "fold${f}: already done, skipping"; continue
  fi
  echo "== fold${f} ${hits[0]}"; date
  singularity exec --nv --env PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages --env TF_FORCE_GPU_ALLOW_GROWTH=true /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif bash -c " chrombpnet contribs_bw       -m ${hits[0]} -r /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/motifscan/peaksets/peaks2k_MEL_IFN.narrowPeak -g /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/data/downloads/hg38.fa -c /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/data/downloads/hg38.chrom.sizes       -op ${op} -pc counts profile"
  date
done
echo "done MEL_IFN"

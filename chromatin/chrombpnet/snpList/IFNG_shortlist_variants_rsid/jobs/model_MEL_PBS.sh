#!/bin/bash
module load conda; source activate fastQTL
set -euo pipefail
shopt -s nullglob
cd /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/variant-scorer/src
for f in 0 1 2 3 4; do
  hits=(/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_PBS_merged_alldonors/biasFactor0.5/fold${f}/chrombpnet_model/models/*chrombpnet_nobias.h5)
  if [[ ${#hits[@]} -ne 1 ]]; then
    echo "ERROR: expected 1 *chrombpnet_nobias.h5 in fold${f}, found ${#hits[@]}" >&2; exit 1
  fi
  model=${hits[0]}
  pre=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/snpList/IFNG_shortlist_variants_rsid/per_fold/MEL_PBS_fold${f}
  echo "== fold${f} ${model}"

  singularity exec --nv --env PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages --env TF_FORCE_GPU_ALLOW_GROWTH=true /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif bash -c "export PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages; python variant_shap.py       -l /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/snpList/IFNG_shortlist_variants_rsid/variants.tsv -g /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/data/downloads/hg38.fa -m ${model} -o ${pre} -s /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/data/downloads/hg38.chrom.sizes       -sc chrombpnet -st counts profile -bs 8"

  singularity exec --nv --env PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages --env TF_FORCE_GPU_ALLOW_GROWTH=true /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif bash -c "export PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages; python variant_scoring.py       -l /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/snpList/IFNG_shortlist_variants_rsid/variants.tsv -g /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/data/downloads/hg38.fa -m ${model} -o ${pre} -s /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/data/downloads/hg38.chrom.sizes       -sc chrombpnet -n 10"
done

# batched ISM: one model load per fold, every variant scored against it
singularity exec --nv --env PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages --env TF_FORCE_GPU_ALLOW_GROWTH=true /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif bash -c "export PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages; python /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/step6b_snpList_ism.py --vs_src /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/variant-scorer/src     --variants /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/snpList/IFNG_shortlist_variants_rsid/variants.tsv --genome /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/data/downloads/hg38.fa     --model_glob '/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_PBS_merged_alldonors/biasFactor0.5/fold{f}/chrombpnet_model/models/*chrombpnet_nobias.h5'     --folds 0,1,2,3,4 --flank 50 --out /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/snpList/IFNG_shortlist_variants_rsid/per_fold/MEL_PBS.ism.tsv"

# Fi-NeMo: all REF+ALT regions of every variant in one fit per head
mh=(/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_PBS_merged_alldonors/biasFactor0.5/fold0/chrombpnet_model/auxiliary/interpret_subsample/*modisco_results_profile_scores.h5)
if [[ ${#mh[@]} -ne 1 ]]; then
  echo "WARNING: expected 1 MoDISco h5 for MEL_PBS, found ${#mh[@]}; skipping Fi-NeMo" >&2
else
  for h in counts profile; do
    fd=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/snpList/IFNG_shortlist_variants_rsid/motifs/finemo/MEL_PBS/${h}
    mkdir -p ${fd}
    h5s=(/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/snpList/IFNG_shortlist_variants_rsid/per_fold/MEL_PBS_fold*.variant_shap.${h}.h5)
    { CUDA_VISIBLE_DEVICES="" /home/shuo.shan-umw/.conda/envs/finemo/bin/finemo extract-regions-chrombpnet-h5 -c "${h5s[@]}" -w 200 -o ${fd}/regions.npz &&
      CUDA_VISIBLE_DEVICES="" /home/shuo.shan-umw/.conda/envs/finemo/bin/finemo call-hits -r ${fd}/regions.npz -m ${mh[0]} -o ${fd} -c 0.0001 -S 10000 &&
      [[ -s ${fd}/hits.tsv ]] &&
      echo "${mh[0]}" > ${fd}/modisco_source.txt; }       || echo "WARNING: Fi-NeMo failed for MEL_PBS ${h}" >&2
  done
fi
echo "done MEL_PBS"

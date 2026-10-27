#!/bin/bash
module load conda; source activate fastQTL
set -euo pipefail
shopt -s nullglob
cd /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/variant-scorer/src
for f in 0 1 2 3 4; do
  # pipeline prefixes model files with <sample>_fold<f>_, so glob for the nobias model
  hits=(/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_FRB_IFN_merged_alldonors/biasFactor0.5/fold${f}/chrombpnet_model/models/*chrombpnet_nobias.h5)
  if [[ ${#hits[@]} -ne 1 ]]; then
    echo "ERROR: expected 1 *chrombpnet_nobias.h5 in fold${f}, found ${#hits[@]}: ${hits[*]}" >&2
    exit 1
  fi
  model=${hits[0]}
  pre=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/posterchild/rs10424564/per_fold/FRB_IFN_fold${f}
  echo "fold${f} model: ${model}"

  # REF/ALT contribution scores (counts + profile heads)
  singularity exec --nv --env PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages --env TF_FORCE_GPU_ALLOW_GROWTH=true /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif bash -c "export PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages; python variant_shap.py       -l /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/posterchild/rs10424564/rs10424564.variant.tsv -g /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/data/downloads/hg38.fa -m ${model} -o ${pre} -s /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/data/downloads/hg38.chrom.sizes       -sc chrombpnet -st counts profile -bs 1"

  # REF/ALT predicted profiles + counts (saved to ${pre}.variant_predictions.h5)
  singularity exec --nv --env PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages --env TF_FORCE_GPU_ALLOW_GROWTH=true /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif bash -c "export PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages; python variant_scoring.py       -l /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/posterchild/rs10424564/rs10424564.variant.tsv -g /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/data/downloads/hg38.fa -m ${model} -o ${pre} -s /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/data/downloads/hg38.chrom.sizes       -sc chrombpnet -n 10"
done

# saturation mutagenesis around the variant, all folds in one python process
singularity exec --nv --env PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages --env TF_FORCE_GPU_ALLOW_GROWTH=true /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif bash -c "export PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages; python /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/step6_ism_snp.py --vs_src /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/variant-scorer/src     --variant /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/posterchild/rs10424564/rs10424564.variant.tsv --genome /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/data/downloads/hg38.fa     --model_glob '/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_FRB_IFN_merged_alldonors/biasFactor0.5/fold{f}/chrombpnet_model/models/*chrombpnet_nobias.h5'     --folds 0,1,2,3,4 --flank 50 --out /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/posterchild/rs10424564/per_fold/FRB_IFN.ism.tsv"

# Fi-NeMo motif hits on fold-averaged REF/ALT contributions (region 0 = REF, 1 = ALT),
# against this model's own MoDISco patterns; failures here should not fail the job
mh=(/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_FRB_IFN_merged_alldonors/biasFactor0.5/fold0/chrombpnet_model/auxiliary/interpret_subsample/*modisco_results_profile_scores.h5)
if [[ ${#mh[@]} -ne 1 ]]; then
  echo "WARNING: expected 1 MoDISco h5 for FRB_IFN, found ${#mh[@]}; skipping Fi-NeMo" >&2
else
  for h in counts profile; do
    fd=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/posterchild/rs10424564/motifs/finemo/FRB_IFN/${h}
    mkdir -p ${fd}
    h5s=(/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/posterchild/rs10424564/per_fold/FRB_IFN_fold*.variant_shap.${h}.h5)
    { CUDA_VISIBLE_DEVICES="" /home/shuo.shan-umw/.conda/envs/finemo/bin/finemo extract-regions-chrombpnet-h5 -c "${h5s[@]}" -w 200 -o ${fd}/regions.npz &&
      CUDA_VISIBLE_DEVICES="" /home/shuo.shan-umw/.conda/envs/finemo/bin/finemo call-hits -r ${fd}/regions.npz -m ${mh[0]} -o ${fd} -c 0.0001 -S 10000 &&
      [[ -s ${fd}/hits.tsv ]] &&
      echo "${mh[0]}" > ${fd}/modisco_source.txt; }       || echo "WARNING: Fi-NeMo failed for FRB_IFN ${h}" >&2
  done
fi
echo "done FRB_IFN"

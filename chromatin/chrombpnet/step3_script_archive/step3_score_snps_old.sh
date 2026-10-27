#!/bin/bash
#BSUB -J chrombpnet_score
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -q long
#BSUB -W 72:00
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/step2_score.%J_%I.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/step2_score.%J_%I.out"

# ============================================================
# DEFINE ONCE HERE — only touch these when changing runs
# ============================================================
BIAS_FACTOR=0.5
FOLD=0
sample=KRT_IFN
prefix=ATAC_${sample}_merged_alldonors

dir="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/${prefix}_biasFactor${BIAS_FACTOR}"
sif="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif"

# ============================================================
# SCRIPT TO SCORE VARIANTS
# ============================================================
rm -f ${dir}/run_inside_singularity_scoreSNPs.sh
cat > ${dir}/run_inside_singularity_scoreSNPs.sh << EOF
#!/bin/bash
dir=${dir}
prefix=${prefix}
BIAS_FACTOR=${BIAS_FACTOR}
FOLD=${FOLD}


chrombpnet score_snps \
    -m chrombpnet_model_fold0/models/chrombpnet_nobias.h5 \
       chrombpnet_model_fold1/models/chrombpnet_nobias.h5 \
       chrombpnet_model_fold2/models/chrombpnet_nobias.h5 \
       chrombpnet_model_fold3/models/chrombpnet_nobias.h5 \
       chrombpnet_model_fold4/models/chrombpnet_nobias.h5 \
    -s variants.bed \
    -g hg38.fa \
    -o variant_scores_output/
echo "Pipeline complete"; date
EOF



echo "Waiting for filesystem sync..."; sleep 10
singularity exec --nv ${sif} bash ${dir}/run_inside_singularity.sh


# output 
# For each variant, you get:
# log_counts_diff: log2(ALT/REF) predicted counts, Quantifies predicted accessibility change
# profile_jsd: JSD between REF and ALT predicted profiles, Detects profile shape changes
# Mean + std across models

#!/bin/bash
base_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_PBS_merged_alldonors/biasFactor0.5
data_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_PBS_merged_alldonors/biasFactor0.5/data
score_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_PBS_merged_alldonors/biasFactor0.5/snp_scores
prefix=ATAC_KRT_PBS_merged_alldonors
snp_file=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap/QTL_sig_snps_chrombpnet_format.bed

echo "================================================"; date
echo "Scoring SNPs for ATAC_KRT_PBS_merged_alldonors"
echo "Total SNPs: $(wc -l < ${snp_file})"
echo "================================================"

# ============================================================
# SCORE SNPS FOR EACH FOLD
# ============================================================
for fold_n in 0 1 2 3 4; do
    model=${base_dir}/fold${fold_n}/chrombpnet_model/models/chrombpnet_nobias.h5

    if [ ! -f "${model}" ]; then
        echo "WARNING: fold ${fold_n} model not found, skipping..."; continue
    fi

    if [ ! -f "${score_dir}/fold${fold_n}_snp_scores.counts_scores.tsv" ]; then
        echo "--- Scoring SNPs for fold ${fold_n} ---"; date
        chrombpnet snp_score             -snps ${snp_file}             -m ${model}             -g ${data_dir}/downloads/hg38.fa             -op ${score_dir}/fold${fold_n}_snp_scores             -bs 64
        echo "--- fold ${fold_n} done ---"; date
    else
        echo "fold ${fold_n} scores already exist, skipping..."
    fi
done

echo "All folds scored"; date
echo "================================================"
echo "Summary:"
echo "  Score files: ${score_dir}/fold*_snp_scores*.tsv"
echo "================================================"; date

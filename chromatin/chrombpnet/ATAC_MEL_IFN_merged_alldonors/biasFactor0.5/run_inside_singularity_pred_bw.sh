#!/bin/bash
base_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5
data_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5/data
pred_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5/pred_bigwigs
prefix=ATAC_MEL_IFN_merged_alldonors

echo "================================================"; date
echo "Predicting bigwigs for ATAC_MEL_IFN_merged_alldonors"
echo "================================================"

# ============================================================
# [1/2] PREDICT BIGWIG FOR EACH FOLD
# ============================================================
for fold_n in 0 1 2 3 4; do
    model=${base_dir}/fold${fold_n}/chrombpnet_model/models/chrombpnet_nobias.h5

    if [ ! -f "${model}" ]; then
        echo "WARNING: fold ${fold_n} model not found at ${model}, skipping..."; continue
    fi

    if [ ! -f "${pred_dir}/fold${fold_n}_chrombpnet_nobias.bw" ]; then
        echo "--- Predicting bigwig for fold ${fold_n} ---"; date
        chrombpnet pred_bw             -cmb ${model}             -g ${data_dir}/downloads/hg38.fa             -c ${data_dir}/downloads/hg38.chrom.sizes             -r ${data_dir}/downloads/peaks_no_blacklist.narrowPeak             -op ${pred_dir}/fold${fold_n}             -bs 64
        echo "--- fold ${fold_n} done ---"; date
    else
        echo "fold ${fold_n} bigwig already exists, skipping..."
    fi
done

echo "All folds predicted"; date

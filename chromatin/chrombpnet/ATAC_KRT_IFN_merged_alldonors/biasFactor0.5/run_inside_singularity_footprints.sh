#!/bin/bash
base_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5
data_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/data
footprint_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/footprints
prefix=ATAC_KRT_IFN_merged_alldonors

echo "================================================"; date
echo "TF footprinting for ATAC_KRT_IFN_merged_alldonors"
echo "TFs to footprint:"; cat ${footprint_dir}/motifs_to_pwm.tsv | cut -f1
echo "================================================"

# ============================================================
# FOOTPRINT FOR EACH FOLD
# ============================================================
for fold_n in 0 1 2 3 4; do
    model=${base_dir}/fold${fold_n}/chrombpnet_model/models/chrombpnet_nobias.h5
    fold_json=${data_dir}/splits/fold_${fold_n}.json

    if [ ! -f "${model}" ]; then
        echo "WARNING: fold ${fold_n} model not found, skipping..."; continue
    fi

    fold_footprint_dir=${footprint_dir}/fold${fold_n}
    mkdir -p ${fold_footprint_dir}

    if [ ! -f "${fold_footprint_dir}/footprints_done.flag" ]; then
        echo "--- Footprinting fold ${fold_n} ---"; date
        chrombpnet footprints             -m ${model}             -g ${data_dir}/downloads/hg38.fa             -r ${data_dir}/downloads/peaks_no_blacklist.narrowPeak             -fl ${fold_json}             -pwm_f ${footprint_dir}/motifs_to_pwm.tsv             -op ${fold_footprint_dir}/footprints             -bs 64
        touch ${fold_footprint_dir}/footprints_done.flag
        echo "--- fold ${fold_n} done ---"; date
    else
        echo "fold ${fold_n} footprints already exist, skipping..."
    fi
done

echo "All folds footprinted"; date
echo "================================================"
echo "Summary:"
echo "  Footprint files: ${footprint_dir}/fold*/footprints*"
echo "================================================"; date

#!/bin/bash
base_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_PBS_merged_alldonors/biasFactor0.6
data_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_PBS_merged_alldonors/biasFactor0.6/data
fold_dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_PBS_merged_alldonors/biasFactor0.6/fold0
prefix=ATAC_KRT_PBS_merged_alldonors
BIAS_FACTOR=0.6
FOLD=0

# fold 0 bias model path — reused for folds 1-4
bias_model_fold0=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_PBS_merged_alldonors/biasFactor0.6/fold0/bias_model/models/ATAC_KRT_PBS_merged_alldonors_fold0_bias.h5

# ============================================================
# FOLD CHROMOSOME DEFINITIONS
# ============================================================
declare -A FOLD_TCR
declare -A FOLD_VCR
FOLD_TCR[0]="chr1 chr3 chr6"
FOLD_VCR[0]="chr8 chr20"
FOLD_TCR[1]="chr2 chr8 chr9 chr16"
FOLD_VCR[1]="chr12 chr17"
FOLD_TCR[2]="chr4 chr11 chr12 chr15 chrY"
FOLD_VCR[2]="chr22 chr7"
FOLD_TCR[3]="chr5 chr10 chr14 chr18 chr20 chr22"
FOLD_VCR[3]="chr6 chr21"
FOLD_TCR[4]="chr7 chr13 chr17 chr19 chr21 chrX"
FOLD_VCR[4]="chr10 chr18"

TCR=${FOLD_TCR[$FOLD]}
VCR=${FOLD_VCR[$FOLD]}

echo "================================================"; date
echo "Fold ${FOLD}: train=${TCR} | val=${VCR}"
echo "================================================"

# ============================================================
echo "--- [1/5] Preparing chromosome splits ---"; date

# write chrom.subset.sizes only if it doesn't exist
if [ ! -f "${data_dir}/downloads/hg38.chrom.subset.sizes" ]; then
    head -n 24 ${data_dir}/downloads/hg38.chrom.sizes         > ${data_dir}/downloads/hg38.chrom.subset.sizes
else
    echo "hg38.chrom.subset.sizes already exists, skipping..."
fi

chrombpnet prep splits   -c ${data_dir}/downloads/hg38.chrom.subset.sizes   -tcr ${TCR}   -vcr ${VCR}   -op ${data_dir}/splits/fold_${FOLD}
echo "--- [1/5] Done ---"; date

# ============================================================
echo "--- [2/5] Cleaning up any partial previous run ---"; date
rm -rf ${fold_dir}/data/nonpeaks_auxiliary/
rm -f  ${fold_dir}/data/nonpeaks_negatives.bed
rm -rf ${fold_dir}/bias_model/logs/
rm -rf ${fold_dir}/bias_model/models/
rm -rf ${fold_dir}/bias_model/auxiliary/
rm -rf ${fold_dir}/bias_model/evaluation/
rm -rf ${fold_dir}/chrombpnet_model/logs/
rm -rf ${fold_dir}/chrombpnet_model/models/
rm -rf ${fold_dir}/chrombpnet_model/auxiliary/
rm -rf ${fold_dir}/chrombpnet_model/evaluation/
echo "--- [2/5] Done ---"; date

# ============================================================
echo "--- [3/5] Generating non-peaks ---"; date
mkdir -p ${fold_dir}/data
chrombpnet prep nonpeaks   -g ${data_dir}/downloads/hg38.fa   -p ${data_dir}/downloads/peaks_no_blacklist.narrowPeak   -c ${data_dir}/downloads/hg38.chrom.sizes   -fl ${data_dir}/splits/fold_${FOLD}.json   -br ${data_dir}/downloads/blacklist.bed.gz   -o ${fold_dir}/data/nonpeaks
echo "--- [3/5] Done ---"; date

# ============================================================
# [4/5] BIAS MODEL — train for fold 0, reuse for folds 1-4
# ============================================================
if [ "${FOLD}" -eq 0 ]; then
    echo "--- [4/5] Training bias model (fold 0) ---"; date
    chrombpnet bias pipeline       -ibam ${data_dir}/downloads/merged.bam       -d "ATAC"       -g ${data_dir}/downloads/hg38.fa       -c ${data_dir}/downloads/hg38.chrom.sizes       -p ${data_dir}/downloads/peaks_no_blacklist.narrowPeak       -n ${fold_dir}/data/nonpeaks_negatives.bed       -fl ${data_dir}/splits/fold_${FOLD}.json       -b ${BIAS_FACTOR}       -o ${fold_dir}/bias_model/       -s 42       -fp ${prefix}_fold${FOLD}

    bias_model=$(ls ${fold_dir}/bias_model/models/*bias.h5 2>/dev/null | head -n 1)
    if [ -z "${bias_model}" ]; then
        echo "ERROR: bias model .h5 not found."; exit 1
    fi
else
    echo "--- [4/5] Reusing bias model from fold 0 ---"; date
    bias_model=${bias_model_fold0}
    if [ ! -f "${bias_model}" ]; then
        echo "ERROR: fold 0 bias model not found at ${bias_model}"
        echo "Make sure fold 0 has completed before running fold ${FOLD}"; exit 1
    fi
fi
echo "--- [4/5] Done — using bias model: ${bias_model} ---"; date

# ============================================================
echo "--- [5/5] Training ChromBPNet model ---"; date
chrombpnet pipeline   -ibam ${data_dir}/downloads/merged.bam   -d "ATAC"   -g ${data_dir}/downloads/hg38.fa   -c ${data_dir}/downloads/hg38.chrom.sizes   -p ${data_dir}/downloads/peaks_no_blacklist.narrowPeak   -n ${fold_dir}/data/nonpeaks_negatives.bed   -fl ${data_dir}/splits/fold_${FOLD}.json   -b ${bias_model}   -o ${fold_dir}/chrombpnet_model/
echo "--- [5/5] Done ---"; date

echo "================================================"
echo "Summary:"
echo "  Fold:             ${FOLD}"
echo "  Bias model:       ${bias_model}"
echo "  ChromBPNet model: ${fold_dir}/chrombpnet_model/"
echo "  Bias factor:      ${BIAS_FACTOR}"
echo "================================================"; date

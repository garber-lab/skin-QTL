#!/bin/bash
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_F110_KRT_IFNG_S1_biasFactor0.5
prefix=ATAC_F110_KRT_IFNG_S1
BIAS_FACTOR=0.5
FOLD=0

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

echo "Fold ${FOLD}: train=${TCR} | val=${VCR}"; date

echo "Preparing chromosome splits..."; date
head -n 24 ${dir}/data/downloads/hg38.chrom.sizes > ${dir}/data/downloads/hg38.chrom.subset.sizes

chrombpnet prep splits   -c ${dir}/data/downloads/hg38.chrom.subset.sizes   -tcr ${TCR}   -vcr ${VCR}   -op ${dir}/data/splits/fold_${FOLD}

echo "Generating non-peaks..."; date
# clean up any partial previous run
rm -rf ${dir}/data/nonpeaks_auxiliary/
rm -f ${dir}/data/nonpeaks_negatives.bed
rm -rf ${dir}/bias_model/logs/
rm -rf ${dir}/bias_model/models/
rm -rf ${dir}/bias_model/auxiliary/
rm -rf ${dir}/bias_model/evaluation/
rm -rf ${dir}/chrombpnet_model/logs/
rm -rf ${dir}/chrombpnet_model/models/
rm -rf ${dir}/chrombpnet_model/auxiliary/
rm -rf ${dir}/chrombpnet_model/evaluation/

chrombpnet prep nonpeaks   -g ${dir}/data/downloads/hg38.fa   -p ${dir}/data/downloads/peaks_no_blacklist.narrowPeak   -c ${dir}/data/downloads/hg38.chrom.sizes   -fl ${dir}/data/splits/fold_${FOLD}.json   -br ${dir}/data/downloads/blacklist.bed.gz   -o ${dir}/data/nonpeaks

echo "Training bias model..."; date
chrombpnet bias pipeline   -ibam ${dir}/data/downloads/merged.bam   -d "ATAC"   -g ${dir}/data/downloads/hg38.fa   -c ${dir}/data/downloads/hg38.chrom.sizes   -p ${dir}/data/downloads/peaks_no_blacklist.narrowPeak   -n ${dir}/data/nonpeaks_negatives.bed   -fl ${dir}/data/splits/fold_${FOLD}.json   -b ${BIAS_FACTOR}   -o ${dir}/bias_model/   -s 42   -fp ${prefix}_fold${FOLD}

bias_model=$(ls ${dir}/bias_model/models/*bias.h5 2>/dev/null | head -n 1)
if [ -z "${bias_model}" ]; then
  echo "ERROR: bias model .h5 not found. Check bias pipeline logs."; exit 1
fi
echo "Using bias model: ${bias_model}"; date

echo "Training ChromBPNet model..."; date
chrombpnet pipeline   -ibam ${dir}/data/downloads/merged.bam   -d "ATAC"   -g ${dir}/data/downloads/hg38.fa   -c ${dir}/data/downloads/hg38.chrom.sizes   -p ${dir}/data/downloads/peaks_no_blacklist.narrowPeak   -n ${dir}/data/nonpeaks_negatives.bed   -fl ${dir}/data/splits/fold_${FOLD}.json   -b ${bias_model}   -o ${dir}/chrombpnet_model/

echo "Pipeline complete"; date

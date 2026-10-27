#!/bin/bash
#BSUB -J chrombpnet_F111_KRT_IFNG_real_vs_pred
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=5000]"
#BSUB -q long
#BSUB -W 72:00
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/step5_F111_KRT_IFNG_real_vs_pred_%J_%I.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/step5_F111_KRT_IFNG_real_vs_pred_%J_%I.out"

# --- define sample paths --- #
module load deeptools

BIAS_FACTOR=0.5
sample=F111_KRT_IFNG
prefix=ATAC_F111_KRT_IFNG_S1
bam_prefix=ATAC_F111_KRT_IFNG_S1

BAM_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/bam_dedupped
BEDPE_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/BEDPE_TAGALIGN_for_chrombpnet
MACS2_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/MACS2_for_chrombpnet
base_dir="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/${prefix}/biasFactor${BIAS_FACTOR}"
data_dir="${base_dir}/data"
pred_dir="${base_dir}/pred_bigwigs"
sif="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif"

compare_dir=${base_dir}/compare_real_vs_predicted
mkdir -p ${compare_dir}
cd ${compare_dir}

## --- actual measurement BAM to bigwig for direct comparison --- #
## counts only the 5' cut site per read, applies the +4/-5 Tn5 shift
#bamCoverage \
#  -b ${BAM_DIR}/${bam_prefix}.nodup.bam \
#  -o ${data_dir}/${prefix}.actual.bw \
#  --normalizeUsing None \
#  --binSize 1 \
#  --Offset 1 1 \
#  -p 8
#
# bigWigInfo ${data_dir}/${prefix}.actual.bw


# --- compare actual to predicted bigwig over the same bed regions --- #
grep -v "^chrM" ${data_dir}/downloads/peaks_no_blacklist.narrowPeak \
  > ${data_dir}/downloads/peaks_no_blacklist_no_chrM.narrowPeak

multiBigwigSummary BED-file \
  --bwfiles \
    ${data_dir}/${prefix}.actual.bw \
    ${pred_dir}/fold0_chrombpnet_nobias.bw \
    ${pred_dir}/fold1_chrombpnet_nobias.bw \
    ${pred_dir}/fold2_chrombpnet_nobias.bw \
    ${pred_dir}/fold3_chrombpnet_nobias.bw \
    ${pred_dir}/fold4_chrombpnet_nobias.bw \
    ${pred_dir}/${prefix}_averaged_chrombpnet_nobias.bw \
  --BED ${data_dir}/downloads/peaks_no_blacklist_no_chrM.narrowPeak \
  --outRawCounts ${compare_dir}/actual_vs_predicted_all_folds_peak_signal.tsv \
  -out ${compare_dir}/actual_vs_predicted_all_folds_peak_signal.npz \
  -p 8

# --- quantify the spearman correlation --- #
plotCorrelation \
  -in ${compare_dir}/actual_vs_predicted_all_folds_peak_signal.npz \
  --corMethod spearman \
  --whatToPlot heatmap \
  --removeOutliers \
  -o ${compare_dir}/actual_vs_predicted_all_folds_spearman_heatmap.pdf \
  --outFileCorMatrix ${compare_dir}/actual_vs_predicted_all_folds_spearman_matrix.tsv

# --- quantify the pearson correlation (log2-transformed for skewed ATAC signal) --- #
plotCorrelation \
  -in ${compare_dir}/actual_vs_predicted_all_folds_peak_signal.npz \
  --corMethod pearson \
  --whatToPlot heatmap \
  --removeOutliers \
  -o ${compare_dir}/actual_vs_predicted_all_folds_pearson_heatmap.pdf \
  --outFileCorMatrix ${compare_dir}/actual_vs_predicted_all_folds_pearson_matrix.tsv


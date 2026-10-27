#!/bin/bash
module load bedtools

workingdir=$1
snp_id=$2
snp_bed="${3:-}"   # empty string if not provided


########## create SNP bed if needed
if [[ ! -f "${snp_bed}" ]]; then
	echo "Warning, SNP bed not provided, generating now"
	bash /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/01_compile_snp_bed.sh ${workingdir} ${snp_id}
	snp_bed=${workingdir}/${snp_id}.bed
fi


########## fetch SNPs that overlap ATACseq peaks
genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
atacPeaks=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/ATACseq_peak_annotation.tsv
atac_bed=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/ATACseq_peak_annotation.bed

cd ${workingdir}

awk 'BEGIN { FS=OFS="\t" } NR == 1 {
    print "snp_chr","snp_start","snp_end","snp",$0
}' "${atacPeaks}" > "${workingdir}/${snp_id}_SNP_ATAC_overlap.txt"

bedtools intersect -a "${snp_bed}" -b "${atac_bed}" -wa -wb >> "${workingdir}/${snp_id}_SNP_ATAC_overlap.txt"

echo "Wrote result to ${workingdir}/${snp_id}_SNP_ATAC_overlap.txt"; date

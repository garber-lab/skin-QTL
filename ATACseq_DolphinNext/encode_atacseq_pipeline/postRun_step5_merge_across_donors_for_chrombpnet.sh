#!/bin/bash
# merge BAMs across all donors, removing chrM
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -W 24:00
#BSUB -q long
#BSUB -J mergeAcrossDonors.KRT
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/log/mergeAcrossDonorsKRT_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/log/mergeAcrossDonorsKRT_%J_%I.err"

module load samtools/1.16.1

# --------------------- #
# Perform cross-donor merging and removal of chrM reads
# --------------------- #
BAM_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/bam_dedupped
cd ${BAM_DIR}

# celltype x condition
sample=KRT_IFN

# ============================================================
# FILTER chrM FROM EACH DONOR BAM
# ============================================================
#for origbamF in ATAC_*${sample}*.nodup.bam; do
#    echo "Filtering chrM from ${origbamF}..."; date
#    samtools view -b -@ 8 \
#        ${origbamF} \
#        chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
#        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
#        chr20 chr21 chr22 chrX chrY \
#        > chrM_removed_${origbamF}
#
#    samtools index chrM_removed_${origbamF}
#
#    echo -n "chrM_removed_${origbamF}: "
#    samtools flagstat chrM_removed_${origbamF} | grep "primary mapped ("
#done
#
# ============================================================
# MERGE ALL DONORS
# ============================================================

# list only .bam files, explicitly exclude .bai
ls chrM_removed_*${sample}*.bam | grep -v "\.bai" > temp_${sample}_files.txt

echo "Files to merge:"; cat temp_${sample}_files.txt

samtools merge -@ 8 -b temp_${sample}_files.txt \
    ATAC_${sample}_merged_alldonors.nodup.bam

samtools index ATAC_${sample}_merged_alldonors.nodup.bam

echo -n "ATAC_${sample}_merged_alldonors.nodup.bam: "
samtools flagstat ATAC_${sample}_merged_alldonors.nodup.bam | grep "primary mapped ("

# ============================================================
# CLEANUP
# ============================================================
rm temp_${sample}_files.txt
rm chrM_removed_*${sample}*.bam   # uncomment once merge is verified

echo "Done"; date

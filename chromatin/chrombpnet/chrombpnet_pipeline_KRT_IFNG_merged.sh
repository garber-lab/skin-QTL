#!/bin/bash
#BSUB -J chrombpnet_merged
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -q gpu
#BSUB -W 720:00
#BSUB -w "done(mergeAcrossDonors.KRT)"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/pipeline.%J_%I.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/pipeline.%J_%I.out"

# ============================================================
# DEFINE ONCE HERE — only touch these when changing runs
# ============================================================
BIAS_FACTOR=0.5
FOLD=0
sample=KRT_IFN
prefix=ATAC_${sample}_merged_alldonors
bam_prefix=${prefix}

BAM_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/bam_dedupped
BEDPE_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/BEDPE_TAGALIGN_for_chrombpnet
MACS2_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/MACS2_for_chrombpnet

dir="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/${prefix}_biasFactor${BIAS_FACTOR}"
tutorial_dir="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial"
sif="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif"

# ============================================================
# SETUP DIRECTORIES
# ============================================================
mkdir -p ${dir}/data/downloads
mkdir -p ${dir}/data/splits
mkdir -p ${dir}/data/output
mkdir -p ${dir}/bias_model
mkdir -p ${dir}/chrombpnet_model
mkdir -p ${dir}/log
mkdir -p ${BEDPE_DIR}
mkdir -p ${MACS2_DIR}

# ============================================================
# COPY REFERENCE FILES
# ============================================================
cp ${tutorial_dir}/data/downloads/blacklist.bed.gz  ${dir}/data/downloads/
cp ${tutorial_dir}/data/downloads/hg38.chrom.sizes  ${dir}/data/downloads/
cp ${tutorial_dir}/data/downloads/hg38.fa           ${dir}/data/downloads/
cp ${tutorial_dir}/data/downloads/hg38.fa.fai       ${dir}/data/downloads/

# ============================================================
# LINK BAM + INDEX
# ============================================================
module load samtools

origbamF="${BAM_DIR}/${bam_prefix}.nodup.bam"

echo "Creating symlink for BAM..."; date
ln -sf ${origbamF} ${dir}/data/downloads/merged.bam

echo "Linking/Indexing BAM..."; date
if [ -f "${origbamF}.bai" ]; then
    ln -sf ${origbamF}.bai ${dir}/data/downloads/merged.bam.bai
else
    samtools index ${dir}/data/downloads/merged.bam
fi
echo "BAM ready"; date

# ============================================================
# ACTIVATE CONDA + LOAD MODULES
# ============================================================
chrsz=${dir}/data/downloads/hg38.chrom.sizes
NPEAKS=300000

source /home/shuo.shan-umw/miniconda3/etc/profile.d/conda.sh
conda activate ATACseq

module load samtools
module load bedtools
module load macs2
module load ucsc_utilities/20240312

merged_bam=${dir}/data/downloads/merged.bam

# ============================================================
# NAME SORT BAM (required for bedtools bamtobed -bedpe)
# ============================================================
if [ ! -f "${dir}/data/downloads/${prefix}.namesorted_for_chrombpnet.tmp.bam" ]; then
    echo "Name-sorting BAM..."; date
    samtools sort -n -@ 8 ${merged_bam} \
        -o ${dir}/data/downloads/${prefix}.namesorted_for_chrombpnet.tmp.bam
else
    echo "Name-sorted BAM already exists, skipping..."; date
fi

# ============================================================
# BAM → BEDPE (NO Tn5 shift — ChromBPNet handles this internally)
# ============================================================
if [ ! -f "${BEDPE_DIR}/${prefix}.bedpe.gz" ]; then
    echo "Creating BEDPE..."; date
    bedtools bamtobed -bedpe -mate1 \
        -i ${dir}/data/downloads/${prefix}.namesorted_for_chrombpnet.tmp.bam | \
        gzip -nc > ${BEDPE_DIR}/${prefix}.bedpe.gz
else
    echo "BEDPE already exists, skipping..."; date
fi

# ============================================================
# BEDPE → tagAlign (NO +4/-5 shift)
# ============================================================
if [ ! -f "${BEDPE_DIR}/${prefix}.tagAlign.gz" ]; then
    echo "Creating tagAlign (unshifted)..."; date
    zcat ${BEDPE_DIR}/${prefix}.bedpe.gz | \
    awk 'BEGIN{OFS="\t"}{
        printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",
        $1,$2,$3,$9,$4,$5,$6,$10
    }' | gzip -nc > ${BEDPE_DIR}/${prefix}.tagAlign.gz
else
    echo "tagAlign already exists, skipping..."; date
fi

# ============================================================
# MACS2 — relaxed peaks for ChromBPNet
# ============================================================
if [ ! -f "${MACS2_DIR}/${prefix}.clipped.narrowPeak.gz" ]; then
    echo "MACS2 callpeak..."; date
    cd ${MACS2_DIR}

    macs2 callpeak \
        -t ${BEDPE_DIR}/${prefix}.tagAlign.gz \
        -f BED \
        -n ${prefix} \
        -g hs \
        -p 0.01 \
        --shift -75 \
        --extsize 150 \
        --nomodel \
        --keep-dup all \
        --call-summits

    echo "Sorting and clipping peaks..."; date
    sort -k8gr,8gr ${MACS2_DIR}/${prefix}_peaks.narrowPeak | \
        awk 'BEGIN{OFS="\t"}{$4="Peak_"NR; print $0}' | \
        head -n ${NPEAKS} | \
        gzip -nc > ${MACS2_DIR}/${prefix}.narrowPeak.gz

    gunzip -c ${MACS2_DIR}/${prefix}.narrowPeak.gz | \
        bedClip /dev/stdin ${chrsz} /dev/stdout | \
        gzip -nc > ${MACS2_DIR}/${prefix}.clipped.narrowPeak.gz

    echo "MACS2 done"; date
else
    echo "Peaks already exist, skipping MACS2..."; date
fi

conda deactivate
echo "Conda ATACseq environment deactivated"; date

# ============================================================
# COPY FINAL PEAK FILE TO CHROMBPNET DIR + CLEANUP
# ============================================================
cp ${MACS2_DIR}/${prefix}.clipped.narrowPeak.gz ${dir}/data/downloads/peaks.clipped.narrowPeak.gz
echo "Peaks ready"; date

# ============================================================
# PREPARE PEAKS — filter blacklist
# ============================================================
peakF="${dir}/data/downloads/peaks.clipped.narrowPeak.gz"

if [ ! -f "${dir}/data/downloads/peaks_no_blacklist.narrowPeak" ]; then
    echo "Processing peaks..."; date

    zcat ${peakF} | sort -k1,1 -k2,2n > ${dir}/data/downloads/orig_peak.bed

    bedtools slop \
        -i ${dir}/data/downloads/blacklist.bed.gz \
        -g ${dir}/data/downloads/hg38.chrom.sizes \
        -b 1057 \
        > ${dir}/data/downloads/temp.bed

    bedtools intersect \
        -v \
        -a ${dir}/data/downloads/orig_peak.bed \
        -b ${dir}/data/downloads/temp.bed \
        > ${dir}/data/downloads/peaks_no_blacklist.narrowPeak

    # cleanup intermediates
    rm -f ${dir}/data/downloads/orig_peak.bed
    rm -f ${dir}/data/downloads/temp.bed

    echo "Done with peak preprocessing"; date
else
    echo "peaks_no_blacklist.narrowPeak already exists, skipping..."; date
fi

# ============================================================
# WRITE INNER SINGULARITY SCRIPT
# ============================================================
rm -f ${dir}/run_inside_singularity.sh
cat > ${dir}/run_inside_singularity.sh << EOF
#!/bin/bash
dir=${dir}
prefix=${prefix}
BIAS_FACTOR=${BIAS_FACTOR}
FOLD=${FOLD}

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

TCR=\${FOLD_TCR[\$FOLD]}
VCR=\${FOLD_VCR[\$FOLD]}

echo "================================================"; date
echo "Fold \${FOLD}: train=\${TCR} | val=\${VCR}"
echo "================================================"

echo "--- [1/6] Preparing chromosome splits ---"; date
head -n 24 \${dir}/data/downloads/hg38.chrom.sizes > \${dir}/data/downloads/hg38.chrom.subset.sizes

chrombpnet prep splits \
  -c \${dir}/data/downloads/hg38.chrom.subset.sizes \
  -tcr \${TCR} \
  -vcr \${VCR} \
  -op \${dir}/data/splits/fold_\${FOLD}
echo "--- [1/6] Done ---"; date

echo "--- [2/6] Cleaning up any partial previous run ---"; date
rm -rf \${dir}/data/nonpeaks_auxiliary/
rm -f  \${dir}/data/nonpeaks_negatives.bed
rm -rf \${dir}/bias_model/logs/
rm -rf \${dir}/bias_model/models/
rm -rf \${dir}/bias_model/auxiliary/
rm -rf \${dir}/bias_model/evaluation/
rm -rf \${dir}/chrombpnet_model/logs/
rm -rf \${dir}/chrombpnet_model/models/
rm -rf \${dir}/chrombpnet_model/auxiliary/
rm -rf \${dir}/chrombpnet_model/evaluation/
echo "--- [2/6] Done ---"; date

echo "--- [3/6] Generating non-peaks ---"; date
chrombpnet prep nonpeaks \
  -g \${dir}/data/downloads/hg38.fa \
  -p \${dir}/data/downloads/peaks_no_blacklist.narrowPeak \
  -c \${dir}/data/downloads/hg38.chrom.sizes \
  -fl \${dir}/data/splits/fold_\${FOLD}.json \
  -br \${dir}/data/downloads/blacklist.bed.gz \
  -o \${dir}/data/nonpeaks
echo "--- [3/6] Done ---"; date

echo "--- [4/6] Training bias model ---"; date
chrombpnet bias pipeline \
  -ibam \${dir}/data/downloads/merged.bam \
  -d "ATAC" \
  -g \${dir}/data/downloads/hg38.fa \
  -c \${dir}/data/downloads/hg38.chrom.sizes \
  -p \${dir}/data/downloads/peaks_no_blacklist.narrowPeak \
  -n \${dir}/data/nonpeaks_negatives.bed \
  -fl \${dir}/data/splits/fold_\${FOLD}.json \
  -b \${BIAS_FACTOR} \
  -o \${dir}/bias_model/ \
  -s 42 \
  -fp \${prefix}_fold\${FOLD}

bias_model=\$(ls \${dir}/bias_model/models/*bias.h5 2>/dev/null | head -n 1)
if [ -z "\${bias_model}" ]; then
  echo "ERROR: bias model .h5 not found. Check bias pipeline logs."; exit 1
fi
echo "--- [4/6] Done — using bias model: \${bias_model} ---"; date

echo "--- [5/6] Training ChromBPNet model ---"; date
chrombpnet pipeline \
  -ibam \${dir}/data/downloads/merged.bam \
  -d "ATAC" \
  -g \${dir}/data/downloads/hg38.fa \
  -c \${dir}/data/downloads/hg38.chrom.sizes \
  -p \${dir}/data/downloads/peaks_no_blacklist.narrowPeak \
  -n \${dir}/data/nonpeaks_negatives.bed \
  -fl \${dir}/data/splits/fold_\${FOLD}.json \
  -b \${bias_model} \
  -o \${dir}/chrombpnet_model/
echo "--- [5/6] Done ---"; date

echo "--- [6/6] Pipeline complete ---"; date
echo "================================================"
echo "Summary:"
echo "  Bias model:       \${bias_model}"
echo "  ChromBPNet model: \${dir}/chrombpnet_model/"
echo "  Fold:             \${FOLD}"
echo "  Bias factor:      \${BIAS_FACTOR}"
echo "================================================"; date
EOF

# ============================================================
# RUN INSIDE SINGULARITY
# ============================================================
echo "Waiting for filesystem sync..."; sleep 10
singularity exec --nv ${sif} bash ${dir}/run_inside_singularity.sh

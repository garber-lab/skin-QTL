#!/bin/bash
#BSUB -J chrombpnet_merged_MEL_IFN_fold4
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4000]"
#BSUB -q gpu
#BSUB -W 72:00
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/pipeline_MEL_IFN_merged_fold4.%J_%I.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/pipeline_MEL_IFN_merged_fold4.%J_%I.out"

# ============================================================
# DEFINE ONCE HERE — only touch these when changing runs
# ============================================================
BIAS_FACTOR=0.5
FOLD=4
sample=MEL_IFN
prefix=ATAC_${sample}_merged_alldonors
bam_prefix=${prefix}

BAM_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/bam_dedupped
BEDPE_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/BEDPE_TAGALIGN_for_chrombpnet
MACS2_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/MACS2_for_chrombpnet

# ============================================================
# DIRECTORY STRUCTURE
# ATAC_KRT_IFN_merged_alldonors/
#   biasFactor0.5/
#     data/          ← shared across all folds
#     fold0/         ← fold-specific
#     fold1/
#     ...
# ============================================================
base_dir="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/${prefix}/biasFactor${BIAS_FACTOR}"
data_dir="${base_dir}/data"
fold_dir="${base_dir}/fold${FOLD}"

tutorial_dir="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial"
sif="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif"

# bias model from fold 0 — reused for folds 1-4
bias_model_fold0="${base_dir}/fold0/bias_model/models/${prefix}_fold0_bias.h5"

# ============================================================
# SETUP DIRECTORIES
# ============================================================
mkdir -p ${data_dir}/downloads
mkdir -p ${data_dir}/splits
mkdir -p ${fold_dir}/bias_model
mkdir -p ${fold_dir}/chrombpnet_model
mkdir -p ${fold_dir}/log
mkdir -p ${BEDPE_DIR}
mkdir -p ${MACS2_DIR}
mkdir -p /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log

# ============================================================
# COPY REFERENCE FILES — skip if already exist (shared across folds)
# ============================================================
echo "--- Checking reference files ---"; date

for f in blacklist.bed.gz hg38.chrom.sizes hg38.fa hg38.fa.fai; do
    if [ ! -f "${data_dir}/downloads/${f}" ]; then
        echo "Copying ${f}..."
        cp ${tutorial_dir}/data/downloads/${f} ${data_dir}/downloads/
    else
        echo "${f} already exists, skipping..."
    fi
done

# ============================================================
# LINK BAM + INDEX — skip if already exists (shared across folds)
# ============================================================
module load samtools

origbamF="${BAM_DIR}/${bam_prefix}.nodup.bam"

if [ ! -f "${data_dir}/downloads/merged.bam" ]; then
    echo "Creating symlink for BAM..."; date
    ln -sf ${origbamF} ${data_dir}/downloads/merged.bam
else
    echo "merged.bam already exists, skipping..."
fi

if [ ! -f "${data_dir}/downloads/merged.bam.bai" ]; then
    echo "Linking/creating BAM index..."; date
    if [ -f "${origbamF}.bai" ]; then
        ln -sf ${origbamF}.bai ${data_dir}/downloads/merged.bam.bai
    else
        samtools index ${data_dir}/downloads/merged.bam
    fi
else
    echo "merged.bam.bai already exists, skipping..."
fi
echo "BAM ready"; date

# ============================================================
# ACTIVATE CONDA + LOAD MODULES
# ============================================================
chrsz=${data_dir}/downloads/hg38.chrom.sizes
NPEAKS=300000

source /home/shuo.shan-umw/miniconda3/etc/profile.d/conda.sh
conda activate ATACseq

module load samtools
module load bedtools
module load macs2
module load ucsc_utilities/20240312

merged_bam=${data_dir}/downloads/merged.bam

# ============================================================
# NAME SORT BAM — skip if already exists (shared across folds)
# ============================================================
if [ ! -f "${data_dir}/downloads/${prefix}.namesorted_for_chrombpnet.tmp.bam" ]; then
    echo "Name-sorting BAM..."; date
    samtools sort -n -@ 8 ${merged_bam} \
        -o ${data_dir}/downloads/${prefix}.namesorted_for_chrombpnet.tmp.bam
else
    echo "Name-sorted BAM already exists, skipping..."; date
fi

# ============================================================
# BAM → BEDPE — skip if already exists (shared across folds)
# ============================================================
if [ ! -f "${BEDPE_DIR}/${prefix}.bedpe.gz" ]; then
    echo "Creating BEDPE..."; date
    bedtools bamtobed -bedpe -mate1 \
        -i ${data_dir}/downloads/${prefix}.namesorted_for_chrombpnet.tmp.bam | \
        gzip -nc > ${BEDPE_DIR}/${prefix}.bedpe.gz
else
    echo "BEDPE already exists, skipping..."; date
fi

# ============================================================
# BEDPE → tagAlign — skip if already exists (shared across folds)
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
# MACS2 — skip if already exists (shared across folds)
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
# COPY PEAK FILE TO SHARED DATA DIR — skip if already exists
# ============================================================
if [ ! -f "${data_dir}/downloads/peaks.clipped.narrowPeak.gz" ]; then
    cp ${MACS2_DIR}/${prefix}.clipped.narrowPeak.gz \
       ${data_dir}/downloads/peaks.clipped.narrowPeak.gz
    echo "Peaks copied to data dir"; date
else
    echo "peaks.clipped.narrowPeak.gz already exists in data dir, skipping..."
fi

# ============================================================
# BLACKLIST FILTER — skip if already exists (shared across folds)
# ============================================================
if [ ! -f "${data_dir}/downloads/peaks_no_blacklist.narrowPeak" ]; then
    echo "Processing peaks — filtering blacklist..."; date

    zcat ${data_dir}/downloads/peaks.clipped.narrowPeak.gz | \
        sort -k1,1 -k2,2n > ${data_dir}/downloads/orig_peak.bed

    bedtools slop \
        -i ${data_dir}/downloads/blacklist.bed.gz \
        -g ${data_dir}/downloads/hg38.chrom.sizes \
        -b 1057 \
        > ${data_dir}/downloads/temp.bed

    bedtools intersect \
        -v \
        -a ${data_dir}/downloads/orig_peak.bed \
        -b ${data_dir}/downloads/temp.bed \
        > ${data_dir}/downloads/peaks_no_blacklist.narrowPeak

    rm -f ${data_dir}/downloads/orig_peak.bed
    rm -f ${data_dir}/downloads/temp.bed

    echo "Done with peak preprocessing"; date
else
    echo "peaks_no_blacklist.narrowPeak already exists, skipping..."; date
fi

# ============================================================
# WRITE INNER SINGULARITY SCRIPT — fold-specific
# ============================================================
rm -f ${fold_dir}/run_inside_singularity.sh
cat > ${fold_dir}/run_inside_singularity.sh << EOF
#!/bin/bash
base_dir=${base_dir}
data_dir=${data_dir}
fold_dir=${fold_dir}
prefix=${prefix}
BIAS_FACTOR=${BIAS_FACTOR}
FOLD=${FOLD}

# fold 0 bias model path — reused for folds 1-4
bias_model_fold0=${bias_model_fold0}

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

# ============================================================
echo "--- [1/5] Preparing chromosome splits ---"; date

# write chrom.subset.sizes only if it doesn't exist
if [ ! -f "\${data_dir}/downloads/hg38.chrom.subset.sizes" ]; then
    head -n 24 \${data_dir}/downloads/hg38.chrom.sizes \
        > \${data_dir}/downloads/hg38.chrom.subset.sizes
else
    echo "hg38.chrom.subset.sizes already exists, skipping..."
fi

chrombpnet prep splits \
  -c \${data_dir}/downloads/hg38.chrom.subset.sizes \
  -tcr \${TCR} \
  -vcr \${VCR} \
  -op \${data_dir}/splits/fold_\${FOLD}
echo "--- [1/5] Done ---"; date

# ============================================================
echo "--- [2/5] Cleaning up any partial previous run ---"; date
rm -rf \${fold_dir}/data/nonpeaks_auxiliary/
rm -f  \${fold_dir}/data/nonpeaks_negatives.bed
rm -rf \${fold_dir}/bias_model/logs/
rm -rf \${fold_dir}/bias_model/models/
rm -rf \${fold_dir}/bias_model/auxiliary/
rm -rf \${fold_dir}/bias_model/evaluation/
rm -rf \${fold_dir}/chrombpnet_model/logs/
rm -rf \${fold_dir}/chrombpnet_model/models/
rm -rf \${fold_dir}/chrombpnet_model/auxiliary/
rm -rf \${fold_dir}/chrombpnet_model/evaluation/
echo "--- [2/5] Done ---"; date

# ============================================================
echo "--- [3/5] Generating non-peaks ---"; date
mkdir -p \${fold_dir}/data
chrombpnet prep nonpeaks \
  -g \${data_dir}/downloads/hg38.fa \
  -p \${data_dir}/downloads/peaks_no_blacklist.narrowPeak \
  -c \${data_dir}/downloads/hg38.chrom.sizes \
  -fl \${data_dir}/splits/fold_\${FOLD}.json \
  -br \${data_dir}/downloads/blacklist.bed.gz \
  -o \${fold_dir}/data/nonpeaks
echo "--- [3/5] Done ---"; date

# ============================================================
# [4/5] BIAS MODEL — train for fold 0, reuse for folds 1-4
# ============================================================
if [ "\${FOLD}" -eq 0 ]; then
    echo "--- [4/5] Training bias model (fold 0) ---"; date
    chrombpnet bias pipeline \
      -ibam \${data_dir}/downloads/merged.bam \
      -d "ATAC" \
      -g \${data_dir}/downloads/hg38.fa \
      -c \${data_dir}/downloads/hg38.chrom.sizes \
      -p \${data_dir}/downloads/peaks_no_blacklist.narrowPeak \
      -n \${fold_dir}/data/nonpeaks_negatives.bed \
      -fl \${data_dir}/splits/fold_\${FOLD}.json \
      -b \${BIAS_FACTOR} \
      -o \${fold_dir}/bias_model/ \
      -s 42 \
      -fp \${prefix}_fold\${FOLD}

    bias_model=\$(ls \${fold_dir}/bias_model/models/*bias.h5 2>/dev/null | head -n 1)
    if [ -z "\${bias_model}" ]; then
        echo "ERROR: bias model .h5 not found."; exit 1
    fi
else
    echo "--- [4/5] Reusing bias model from fold 0 ---"; date
    bias_model=\${bias_model_fold0}
    if [ ! -f "\${bias_model}" ]; then
        echo "ERROR: fold 0 bias model not found at \${bias_model}"
        echo "Make sure fold 0 has completed before running fold \${FOLD}"; exit 1
    fi
fi
echo "--- [4/5] Done — using bias model: \${bias_model} ---"; date

# ============================================================
echo "--- [5/5] Training ChromBPNet model ---"; date
chrombpnet pipeline \
  -ibam \${data_dir}/downloads/merged.bam \
  -d "ATAC" \
  -g \${data_dir}/downloads/hg38.fa \
  -c \${data_dir}/downloads/hg38.chrom.sizes \
  -p \${data_dir}/downloads/peaks_no_blacklist.narrowPeak \
  -n \${fold_dir}/data/nonpeaks_negatives.bed \
  -fl \${data_dir}/splits/fold_\${FOLD}.json \
  -b \${bias_model} \
  -o \${fold_dir}/chrombpnet_model/
echo "--- [5/5] Done ---"; date

echo "================================================"
echo "Summary:"
echo "  Fold:             \${FOLD}"
echo "  Bias model:       \${bias_model}"
echo "  ChromBPNet model: \${fold_dir}/chrombpnet_model/"
echo "  Bias factor:      \${BIAS_FACTOR}"
echo "================================================"; date
EOF

# ============================================================
# RUN INSIDE SINGULARITY
# ============================================================
echo "Waiting for filesystem sync..."; sleep 10
singularity exec --nv ${sif} bash ${fold_dir}/run_inside_singularity.sh

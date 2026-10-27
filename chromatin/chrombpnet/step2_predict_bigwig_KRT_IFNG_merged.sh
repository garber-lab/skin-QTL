#!/bin/bash
#BSUB -J chrombpnet_pred_bw_KRT_IFNG
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -q gpu
#BSUB -W 72:00
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/step2_pred_bw_KRT_IFNG.%J_%I.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/log/step2_pred_bw_KRT_IFNG.%J_%I.out"

# ============================================================
# DEFINE ONCE HERE — only touch these when changing runs
# ============================================================
BIAS_FACTOR=0.5
sample=KRT_IFN
prefix=ATAC_${sample}_merged_alldonors

base_dir="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/${prefix}/biasFactor${BIAS_FACTOR}"
data_dir="${base_dir}/data"
pred_dir="${base_dir}/pred_bigwigs"
sif="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif"

mkdir -p ${pred_dir}

# ============================================================
# WRITE INNER SINGULARITY SCRIPT
# ============================================================
rm -f ${base_dir}/run_inside_singularity_pred_bw.sh
cat > ${base_dir}/run_inside_singularity_pred_bw.sh << EOF
#!/bin/bash
base_dir=${base_dir}
data_dir=${data_dir}
pred_dir=${pred_dir}
prefix=${prefix}

echo "================================================"; date
echo "Predicting bigwigs for ${prefix}"
echo "================================================"

# ============================================================
# [1/2] PREDICT BIGWIG FOR EACH FOLD
# ============================================================
for fold_n in 0 1 2 3 4; do
    model=\${base_dir}/fold\${fold_n}/chrombpnet_model/models/chrombpnet_nobias.h5

    if [ ! -f "\${model}" ]; then
        echo "WARNING: fold \${fold_n} model not found at \${model}, skipping..."; continue
    fi

    if [ ! -f "\${pred_dir}/fold\${fold_n}_chrombpnet_nobias.bw" ]; then
        echo "--- Predicting bigwig for fold \${fold_n} ---"; date
        chrombpnet pred_bw \
            -cmb \${model} \
            -g \${data_dir}/downloads/hg38.fa \
            -c \${data_dir}/downloads/hg38.chrom.sizes \
            -r \${data_dir}/downloads/peaks_no_blacklist.narrowPeak \
            -op \${pred_dir}/fold\${fold_n} \
            -bs 64
        echo "--- fold \${fold_n} done ---"; date
    else
        echo "fold \${fold_n} bigwig already exists, skipping..."
    fi
done

echo "All folds predicted"; date
EOF

# ============================================================
# RUN INSIDE SINGULARITY
# ============================================================
echo "Waiting for filesystem sync..."; date; sleep 10
singularity exec --nv ${sif} bash ${base_dir}/run_inside_singularity_pred_bw.sh

# ============================================================
# AVERAGE BIGWIGS ACROSS FOLDS — outside singularity with deeptools
# ============================================================
echo "--- Averaging bigwigs across folds ---"; date

module load deeptools

# collect bigwig files
bw_files=""
for fold_n in 0 1 2 3 4; do
    bw=${pred_dir}/fold${fold_n}_chrombpnet_nobias.bw
    if [ -f "${bw}" ]; then
        bw_files="${bw_files} ${bw}"
    else
        echo "WARNING: ${bw} not found, skipping from average"
    fi
done

echo "Averaging: ${bw_files}"; date

bigwigAverage \
    -b ${bw_files} \
    -o ${pred_dir}/${prefix}_averaged_chrombpnet_nobias.bw \
    --binSize 1 \
    -p 8 \
    --skipNonCoveredRegions

echo "================================================"
echo "Done!"
echo "  Per-fold bigwigs: ${pred_dir}/fold*_chrombpnet_nobias.bw"
echo "  Averaged bigwig:  ${pred_dir}/${prefix}_averaged_chrombpnet_nobias.bw"
echo "================================================"; date

#!/bin/bash
#BSUB -J ABC_FRB_IFNG
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=32000]"
#BSUB -W 12:00
#BSUB -q long
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/ABC/logs/ABC_FRB_IFNG_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/ABC/logs/ABC_FRB_IFNG_%J.err"
# ============================================================
#  ABC enhancer-gene prediction — context FRB_IFNG
#  Run:  bsub < run_ABC_FRB_IFNG.sh
# ============================================================

# ---------- CONFIG (edit paths here) ----------
CONTEXT=FRB_IFNG
CELLTYPE=FRB
CONDITION=IFNG

# Inputs
ATAC_BAM="/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/bam_dedupped/ATAC_FRB_IFNG_merged_alldonors.nodup.bam"   # merged-donor ATAC bam for this context
H3K27AC_DONOR_BAMS="/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DolphinNext/FRB/report5357/chip/bam/ChIP_F25_FRB_IFN.bam /pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DolphinNext/FRB/report5357/chip/bam/ChIP_F49_FRB_IFN.bam /pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DolphinNext/FRB/report5357/chip/bam/ChIP_F55_FRB_IFN.bam "                                                  # 3 donor H3K27ac bams (merged below)
H3K27AC_BAM="/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/ABC/FRB_IFNG/h3k27ac_FRB_IFNG_merged.bam"                         # produced by the merge step
EXPRESSED_GENES="/pi/manuel.garber-umw/human/skin/eQTLs/ABC/expressed_genes/expressed_FRB_IFNG.txt"                                 # TSV: gene_symbol <TAB> value (see README); or comment out --expression_table below

# Shared references
ABC_REPO="/pi/manuel.garber-umw/human/skin/eQTLs/tools/ABC-Enhancer-Gene-Prediction"
ABC_ENV="abcmodel"
CHROM_SIZES="/share/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.chrom.sizes"
BLACKLIST="/pi/manuel.garber-umw/human/skin/eQTLs/literature/blacklist/hg38.blacklist.bed"
GENES_BED="/pi/manuel.garber-umw/human/skin/eQTLs/annotation/gencode_v34_gene_bounds.bed"
UBIQ_GENES="${ABC_REPO}/reference/UbiquitouslyExpressedGenesHG19.txt"
QNORM="${ABC_REPO}/src/EnhancersQNormRef.K562.txt"

# Contact model:  powerlaw  (no Hi-C, fast)   OR   hic  (set HIC_DIR)
CONTACT_MODE="powerlaw"
HIC_DIR="/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/ABC/hic/FRB"        # only used if CONTACT_MODE=hic ; ABC juicer-dumped dir, 5kb KR-normalized
HIC_RESOLUTION=5000

# Output
OUTDIR="/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/ABC/FRB_IFNG"
# ---------- END CONFIG ----------

set -euo pipefail
mkdir -p "${OUTDIR}/Peaks" "${OUTDIR}/Neighborhoods" "${OUTDIR}/Predictions" "$(dirname "$OUTDIR")/logs"

# conda + set -u interplay (binutils activation hooks use unbound vars)
set +u; source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${ABC_ENV}"; set -u
module load samtools/1.16.1 bedtools/2.30.0 macs2/2.2.7.1 2>/dev/null || true

log(){ echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ---------------------------------------------------------
# Step 0: merge + index the 3 donor H3K27ac bams for this context (skip if exists)
# ---------------------------------------------------------
if [ ! -s "${H3K27AC_BAM}" ]; then
  log "Merging H3K27ac donor bams -> ${H3K27AC_BAM}"
  samtools merge -f -@ 4 "${H3K27AC_BAM}" ${H3K27AC_DONOR_BAMS}
  samtools index "${H3K27AC_BAM}"
else log "H3K27ac merged bam exists, skipping"; fi
[ -s "${ATAC_BAM}.bai" ] || samtools index "${ATAC_BAM}"

# ---------------------------------------------------------
# Step 1: MACS2 peaks -> sorted narrowPeak
# ---------------------------------------------------------
PK="${OUTDIR}/Peaks/${CONTEXT}"
if [ ! -s "${PK}_peaks.narrowPeak.sorted" ]; then
  log "MACS2 callpeak"
  macs2 callpeak -t "${ATAC_BAM}" -n "${CONTEXT}" -f BAM -g hs -p 0.1 --call-summits --outdir "${OUTDIR}/Peaks"
  bedtools sort -faidx "${CHROM_SIZES}" -i "${PK}_peaks.narrowPeak" > "${PK}_peaks.narrowPeak.sorted"
else log "peaks exist, skipping macs2"; fi

# ---------------------------------------------------------
# Step 2: candidate enhancer regions (top 150k, +/-250bp of summit)
# ---------------------------------------------------------
CAND="${PK}_peaks.narrowPeak.sorted.candidateRegions.bed"
if [ ! -s "${CAND}" ]; then
  log "makeCandidateRegions"
  python "${ABC_REPO}/src/makeCandidateRegions.py" \
    --narrowPeak "${PK}_peaks.narrowPeak.sorted" \
    --bam "${ATAC_BAM}" \
    --outDir "${OUTDIR}/Peaks" \
    --chrom_sizes "${CHROM_SIZES}" \
    --regions_blacklist "${BLACKLIST}" \
    --regions_whitelist "${GENES_BED}" \
    --peakExtendFromSummit 250 \
    --nStrongestPeaks 150000
else log "candidate regions exist, skipping"; fi

# ---------------------------------------------------------
# Step 3: neighborhoods (quantify ATAC + H3K27ac; qnorm for cross-context comparability)
# ---------------------------------------------------------
if [ ! -s "${OUTDIR}/Neighborhoods/EnhancerList.txt" ]; then
  log "run.neighborhoods"
  EXPR_FLAG=""; [ -s "${EXPRESSED_GENES}" ] && EXPR_FLAG="--expression_table ${EXPRESSED_GENES}"
  python "${ABC_REPO}/src/run.neighborhoods.py" \
    --candidate_enhancer_regions "${CAND}" \
    --genes "${GENES_BED}" \
    --H3K27ac "${H3K27AC_BAM}" \
    --DHS "${ATAC_BAM}" \
    ${EXPR_FLAG} \
    --chrom_sizes "${CHROM_SIZES}" \
    --ubiquitously_expressed_genes "${UBIQ_GENES}" \
    --qnorm "${QNORM}" \
    --cellType "${CONTEXT}" \
    --outdir "${OUTDIR}/Neighborhoods"
else log "neighborhoods exist, skipping"; fi

# ---------------------------------------------------------
# Step 4: predict ABC scores (power-law OR Hi-C), threshold 0.02
# ---------------------------------------------------------
log "predict (CONTACT_MODE=${CONTACT_MODE})"
if [ "${CONTACT_MODE}" = "hic" ]; then
  python "${ABC_REPO}/src/predict.py" \
    --enhancers "${OUTDIR}/Neighborhoods/EnhancerList.txt" \
    --genes "${OUTDIR}/Neighborhoods/GeneList.txt" \
    --HiCdir "${HIC_DIR}" --hic_resolution "${HIC_RESOLUTION}" --scale_hic_using_powerlaw \
    --threshold 0.02 --cellType "${CONTEXT}" --outdir "${OUTDIR}/Predictions" --make_all_putative
else
  python "${ABC_REPO}/src/predict.py" \
    --enhancers "${OUTDIR}/Neighborhoods/EnhancerList.txt" \
    --genes "${OUTDIR}/Neighborhoods/GeneList.txt" \
    --score_column powerlaw.Score \
    --threshold 0.02 --cellType "${CONTEXT}" --outdir "${OUTDIR}/Predictions" --make_all_putative
fi

log "DONE ${CONTEXT}. Predictions: ${OUTDIR}/Predictions/EnhancerPredictions.txt"

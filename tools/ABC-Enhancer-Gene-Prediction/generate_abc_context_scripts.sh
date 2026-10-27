#!/bin/bash
# ============================================================
#  Generate one self-contained ABC BSUB job script per context
#  (6 contexts = 3 celltypes x 2 conditions). Each output script
#  bakes the context into its name (run_ABC_<CONTEXT>.sh) and has
#  ALL file-path assumptions in a CONFIG block at the top so you
#  can edit paths quickly and submit them in parallel:
#      for f in scripts/run_ABC_*.sh; do bsub < "$f"; done
#
#  Contact defaults to POWER-LAW (no Hi-C download needed) so you
#  can get results fast; switch to Hi-C by editing CONTACT_MODE.
# ============================================================
set -euo pipefail
OUTDIR="$(cd "$(dirname "$0")" && pwd)/scripts"
mkdir -p "$OUTDIR"

# ---- SHARED paths (edit here once; injected into every script) ----
ABC_REPO='/pi/manuel.garber-umw/human/skin/eQTLs/tools/ABC-Enhancer-Gene-Prediction'   # git clone -b v0.2.2 ...
ABC_ENV='abcmodel'                                                                       # conda env with ABC deps (macs2, samtools, pandas, pyBigWig, juicer for hic)
CHROM_SIZES='/share/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.chrom.sizes'
BLACKLIST='/pi/manuel.garber-umw/human/skin/eQTLs/literature/blacklist/hg38.blacklist.bed'
GENES_BED='/pi/manuel.garber-umw/human/skin/eQTLs/annotation/gencode_v34_gene_bounds.bed'   # bed6: chr,start,end,gene_symbol,score,strand (gene bounds; ABC derives TSS)
UBIQ_GENES='${ABC_REPO}/reference/UbiquitouslyExpressedGenesHG19.txt'                       # gene-symbol list, genome-agnostic
QNORM='${ABC_REPO}/src/EnhancersQNormRef.K562.txt'                                          # cross-celltype comparability (IMPORTANT: you compare across contexts)
OUT_ROOT='/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/ABC'
ATAC_DIR='/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/bam_dedupped'
CHIP_ROOT='/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DolphinNext'
EXPR_DIR='/pi/manuel.garber-umw/human/skin/eQTLs/ABC/expressed_genes'                       # per-context expressed-gene tables (see README for format)
CONTACT_MODE='powerlaw'   # 'powerlaw' (no Hi-C) | 'hic' (set HIC_DIR per script)

# ---- per-context table: CTX  CELLTYPE  COND  ATAC_COND  CHIP_REPORT  CHIP_COND  CHIP_PREFIX ----
# CHIP file naming differs: KRT uses  ChIP_F##K_<cond>.bam ; FRB/MEL use ChIP_F##_<CT>_<cond>.bam
# ChIP condition token is IFN (not IFNG) / PBS.
read -r -d '' TABLE <<'TAB' || true
FRB_PBS   FRB PBS  PBS  report5357 PBS  FRB
FRB_IFNG  FRB IFNG IFNG report5357 IFN  FRB
KRT_PBS   KRT PBS  PBS  report5355 PBS  KRT
KRT_IFNG  KRT IFNG IFNG report5355 IFN  KRT
MEL_PBS   MEL PBS  PBS  report5356 PBS  MEL
MEL_IFNG  MEL IFNG IFNG report5356 IFN  MEL
TAB

DONORS="F25 F49 F55"

while read -r CTX CELLTYPE COND ATAC_COND CHIP_REPORT CHIP_COND CHIP_PREFIX; do
  [ -z "${CTX:-}" ] && continue
  CHIP_BAM_DIR="${CHIP_ROOT}/${CELLTYPE}/${CHIP_REPORT}/chip/bam"
  # build the 3 donor H3K27ac bam paths for this context
  CHIP_LIST=""
  for D in $DONORS; do
    if [ "$CELLTYPE" = "KRT" ]; then f="${CHIP_BAM_DIR}/ChIP_${D}K_${CHIP_COND}.bam"
    else f="${CHIP_BAM_DIR}/ChIP_${D}_${CELLTYPE}_${CHIP_COND}.bam"; fi
    CHIP_LIST="${CHIP_LIST}${f} "
  done
  OUT="${OUTDIR}/run_ABC_${CTX}.sh"

  # ---- header + CONFIG (expanded now) ----
  cat > "$OUT" <<EOF
#!/bin/bash
#BSUB -J ABC_${CTX}
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=32000]"
#BSUB -W 12:00
#BSUB -q long
#BSUB -o "${OUT_ROOT}/logs/ABC_${CTX}_%J.out"
#BSUB -e "${OUT_ROOT}/logs/ABC_${CTX}_%J.err"
# ============================================================
#  ABC enhancer-gene prediction — context ${CTX}
#  Run:  bsub < run_ABC_${CTX}.sh
# ============================================================

# ---------- CONFIG (edit paths here) ----------
CONTEXT=${CTX}
CELLTYPE=${CELLTYPE}
CONDITION=${COND}

# Inputs
ATAC_BAM="${ATAC_DIR}/ATAC_${CELLTYPE}_${ATAC_COND}_merged_alldonors.nodup.bam"   # merged-donor ATAC bam for this context
H3K27AC_DONOR_BAMS="${CHIP_LIST}"                                                  # 3 donor H3K27ac bams (merged below)
H3K27AC_BAM="${OUT_ROOT}/${CTX}/h3k27ac_${CTX}_merged.bam"                         # produced by the merge step
EXPRESSED_GENES="${EXPR_DIR}/expressed_${CTX}.txt"                                 # TSV: gene_symbol <TAB> value (see README); or comment out --expression_table below

# Shared references
ABC_REPO="${ABC_REPO}"
ABC_ENV="${ABC_ENV}"
CHROM_SIZES="${CHROM_SIZES}"
BLACKLIST="${BLACKLIST}"
GENES_BED="${GENES_BED}"
UBIQ_GENES="${UBIQ_GENES}"
QNORM="${QNORM}"

# Contact model:  powerlaw  (no Hi-C, fast)   OR   hic  (set HIC_DIR)
CONTACT_MODE="${CONTACT_MODE}"
HIC_DIR="${OUT_ROOT}/hic/${CELLTYPE}"        # only used if CONTACT_MODE=hic ; ABC juicer-dumped dir, 5kb KR-normalized
HIC_RESOLUTION=5000

# Output
OUTDIR="${OUT_ROOT}/${CTX}"
# ---------- END CONFIG ----------
EOF

  # ---- body (literal; uses the shell vars above) ----
  cat >> "$OUT" <<'EOF'

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
EOF
  chmod +x "$OUT"
  echo "wrote $OUT"
done <<< "$TABLE"

echo "Done. Submit with:  for f in ${OUTDIR}/run_ABC_*.sh; do bsub < \"\$f\"; done"

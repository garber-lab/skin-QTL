#!/bin/bash
# ==============================================================================
# Phase 0B -- phyloP 100-way conservation for the Phase 0 variant set
#
# Approach: download the UCSC hg38 phyloP100way bigWig and the UCSC
#           bigWigAverageOverBed binary, then score the variant BED directly.
#           No conversion to bedGraph is needed.
#
# Submit:   bsub < 01_phylop.lsf
# Input:    ${BED}   (variants_phylop.bed, headerless BED4, name = variant_id)
# Output:   ${OUTDIR}/phylop_per_variant.tab
#           columns: name, size, covered, sum, mean0, mean
#           -> use `mean` as the variant's phyloP; `covered` = 0 means the
#              variant has no phyloP data and should be NA, not 0.
# ==============================================================================

#BSUB -J phylop
#BSUB -q long
#BSUB -n 2
#BSUB -R rusage[mem=8000]
#BSUB -W 12:00
#BSUB -o phylop.%J.out
#BSUB -e phylop.%J.err

# ---- EDIT THESE --------------------------------------------------------------
BED=/pi/manuel.garber-umw/human/skin/eQTLs/literature/phyloP/variants_phylop.bed
OUTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/literature/phyloP

BW_URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw
TOOL_URL=https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed
# ------------------------------------------------------------------------------

set -e
cd "${OUTDIR}"

# ~9 GB, downloaded once. Skip if already present.
wget -c "${BW_URL}"

module load ucsc_utilities/20240312
bigWigAverageOverBed \
  hg38.phyloP100way.bw \
  "${BED}" \
  "${OUTDIR}/phylop_per_variant.tab"

echo "done: ${OUTDIR}/phylop_per_variant.tab"
wc -l "${OUTDIR}/phylop_per_variant.tab"


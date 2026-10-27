#!/bin/bash
# ==============================================================================
# Phase 0B -- CADD v1.7 scores for the Phase 0 variant set
#
# CADD is NOT installed on the cluster. Rather than install CADD-scripts (which
# would re-run the whole annotation pipeline), this looks the variants up in
# CADD's precomputed tabix-indexed files. Far simpler and much faster.
#
# Two files are used:
#   whole_genome_SNVs.tsv.gz          ~80 GB, covers every possible SNV
#   gnomad.genomes.r4.0.indel.tsv.gz  ~2 GB, covers gnomAD-observed indels only
# Indels absent from the second file will simply have no CADD score. That is
# expected and should be recorded as NA, not 0.
#
# Submit:   bsub < 03_cadd.lsf
# Input:    ${REGIONS} (variants_regions.bed, sorted BED3)
# Output:   ${OUTDIR}/cadd_snv_hits.tsv    chrom pos ref alt raw phred
#           ${OUTDIR}/cadd_indel_hits.tsv  same columns
#
# ALTERNATIVE, if the 80 GB download is unwelcome: the CADD web service at
# https://cadd.bihealth.org/score accepts a VCF upload of this size and returns
# the same scores. It is much lighter, but it does send your variant
# coordinates to an external server -- your call, not something I will do.
# ==============================================================================

#BSUB -J cadd
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=8000]
#BSUB -W 24:00
#BSUB -q long
#BSUB -o cadd.%J.out
#BSUB -e cadd.%J.err

# ---- EDIT THESE --------------------------------------------------------------
WORKDIR=/pi/manuel.garber-umw/human/skin/eQTLs/literature/CADD
REGIONS=${WORKDIR}/variants_regions.bed
OUTDIR=${WORKDIR}
REFDIR=${WORKDIR}/refdata
CONDA_ENV=htslib

CADD_BASE=https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38
SNV_FILE=whole_genome_SNVs.tsv.gz
INDEL_FILE=gnomad.genomes.r4.0.indel.tsv.gz
# ------------------------------------------------------------------------------

set -e
mkdir -p "${OUTDIR}" "${REFDIR}"

source "$(conda info --base)/etc/profile.d/conda.sh"
if ! conda env list | grep -q "^${CONDA_ENV} "; then
  conda create -y -n "${CONDA_ENV}" -c conda-forge -c bioconda htslib
fi
conda activate "${CONDA_ENV}"

cd "${REFDIR}"

# ---- one-time downloads ------------------------------------------------------
for f in "${SNV_FILE}" "${SNV_FILE}.tbi" "${INDEL_FILE}" "${INDEL_FILE}.tbi"; do
  if [ ! -f "$f" ]; then
    wget -c "${CADD_BASE}/$f"
  fi
done

# ---- lookup ------------------------------------------------------------------
# CADD files use chromosome names WITHOUT the "chr" prefix, so strip it.
sed 's/^chr//' "${REGIONS}" > "${OUTDIR}/regions_nochr.bed"

echo -e "chrom\tpos\tref\talt\tCADD_raw\tCADD_PHRED" > "${OUTDIR}/cadd_snv_hits.tsv"
tabix -R "${OUTDIR}/regions_nochr.bed" "${SNV_FILE}" \
  >> "${OUTDIR}/cadd_snv_hits.tsv"

echo -e "chrom\tpos\tref\talt\tCADD_raw\tCADD_PHRED" > "${OUTDIR}/cadd_indel_hits.tsv"
tabix -R "${OUTDIR}/regions_nochr.bed" "${INDEL_FILE}" \
  >> "${OUTDIR}/cadd_indel_hits.tsv"

echo "done"
wc -l "${OUTDIR}/cadd_snv_hits.tsv" "${OUTDIR}/cadd_indel_hits.tsv"


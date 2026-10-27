 #!/bin/bash
# ==============================================================================
# Phase 0B -- Ensembl VEP v110 consequences for the Phase 0 variant set
#
# VEP and its cache are NOT installed on the cluster, so this script installs
# both via conda on first run. The cache download is ~25 GB and takes a while;
# it is done once and reused.
#
# Submit:   bsub < 02_vep.lsf
# Input:    ${VCF}  (variants_hg38.vcf, 77543 records)
# Output:   ${OUTDIR}/vep_annotations.tsv   tab-delimited, ALL transcripts
#
# NOTE ON --pick: this script deliberately does NOT use --pick, so every
# transcript consequence is reported. Phase 0 needs to know whether ANY
# transcript gives a coding / protein-altering / splicing consequence, and
# --pick would hide a coding consequence behind a chosen non-coding one.
# The trade-off is a larger file (roughly 1-4 M rows). If you want one row per
# variant instead, add --pick_allele_gene to VEP_OPTS below.
# ==============================================================================

#BSUB -J vep110
#BSUB -q long
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=8000]
#BSUB -W 24:00
#BSUB -o vep.%J.out
#BSUB -e vep.%J.err

# ---- EDIT THESE --------------------------------------------------------------
WORKDIR=/pi/manuel.garber-umw/human/skin/eQTLs/literature/VEP
VCF=${WORKDIR}/variants_hg38.vcf
OUTDIR=${WORKDIR}
CACHE_DIR=${WORKDIR}/vep_cache
CONDA_ENV=vep110

VEP_VERSION=110
ASSEMBLY=GRCh38
NFORK=4
# ------------------------------------------------------------------------------

set -e
mkdir -p "${OUTDIR}" "${CACHE_DIR}"

# ---- one-time install --------------------------------------------------------
source "$(conda info --base)/etc/profile.d/conda.sh"

if ! conda env list | grep -q "^${CONDA_ENV} "; then
  conda create -y -n "${CONDA_ENV}" -c conda-forge -c bioconda \
    ensembl-vep=${VEP_VERSION}
fi
conda activate "${CONDA_ENV}"

# ---- one-time cache download (~25 GB) ----------------------------------------
if [ ! -d "${CACHE_DIR}/homo_sapiens/${VEP_VERSION}_${ASSEMBLY}" ]; then
  vep_install -a cf -s homo_sapiens -y ${ASSEMBLY} \
    -c "${CACHE_DIR}" --CONVERT --NO_UPDATE
fi

# ---- run ---------------------------------------------------------------------
vep \
  --input_file  "${VCF}" \
  --output_file "${OUTDIR}/vep_annotations.tsv" \
  --species homo_sapiens \
  --assembly ${ASSEMBLY} \
  --cache --dir_cache "${CACHE_DIR}" --offline \
  --tab \
  --fields "Uploaded_variation,Location,Allele,Gene,SYMBOL,Feature,Feature_type,Consequence,IMPACT,BIOTYPE,CANONICAL,DISTANCE,STRAND" \
  --symbol --canonical --biotype --distance 5000 \
  --fork ${NFORK} \
  --no_stats --force_overwrite

echo "done: ${OUTDIR}/vep_annotations.tsv"
grep -vc '^#' "${OUTDIR}/vep_annotations.tsv"


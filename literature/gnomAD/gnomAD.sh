#!/bin/bash
# ==============================================================================
# Phase 0B -- gnomAD v4.0 allele frequencies for the Phase 0 variant set
#
# The full gnomAD v4.0 genomes sites VCFs are terabytes, so this does NOT
# download them. bcftools reads the public Google bucket over HTTPS and pulls
# only the byte ranges covering our 77,543 positions, one chromosome at a time.
# Nothing large is stored on the cluster.
#
# Submit:   bsub < 04_gnomad.lsf
# Input:    ${REGIONS} (variants_regions.bed, sorted BED3, "chr"-prefixed)
# Output:   ${OUTDIR}/gnomad_af.tsv
#           columns: chrom pos ref alt AF AC AN AF_grpmax AF_nfe AF_afr AF_eas AF_sas AF_amr
#
# REQUIREMENT: the compute nodes must be able to reach the internet over HTTPS.
# If outbound traffic is blocked on compute nodes, run this on a login or data
# transfer node instead:  bash 04_gnomad.lsf
# ==============================================================================

#BSUB -J gnomad
#BSUB -q long
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=4000]
#BSUB -W 24:00
#BSUB -o gnomad.%J.out
#BSUB -e gnomad.%J.err

# ---- EDIT THESE --------------------------------------------------------------
WORKDIR=/pi/manuel.garber-umw/human/skin/eQTLs/literature/gnomAD
REGIONS=/pi/manuel.garber-umw/human/skin/eQTLs/literature/CADD/variants_regions.bed
OUTDIR=${WORKDIR}
BCFTOOLS_MODULE=bcftools        # `module avail bcftools` to see exact name/version

GNOMAD_BASE=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/genomes
GNOMAD_PREFIX=gnomad.genomes.v4.0.sites
CHROMS=$(seq 1 22)
# ------------------------------------------------------------------------------

set -e
mkdir -p "${OUTDIR}/per_chrom"

module load ${BCFTOOLS_MODULE}

# This job reads gnomAD straight off an HTTPS URL, so bcftools must have been
# built against libcurl. Most module builds are, but not all -- and without it
# the failure is a confusing "protocol not supported" partway in. Check once,
# cheaply, before looping over 22 chromosomes.
bcftools view -h "${GNOMAD_BASE}/${GNOMAD_PREFIX}.chr22.vcf.bgz" > /dev/null

# gnomAD v4 on GRCh38 uses "chr"-prefixed contig names, same as our regions
# file, so no renaming is needed here (unlike CADD, which has no prefix).

FMT='%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%INFO/AF_grpmax\t%INFO/AF_nfe\t%INFO/AF_afr\t%INFO/AF_eas\t%INFO/AF_sas\t%INFO/AF_amr\n'

for c in ${CHROMS}; do
  echo "=== chr${c} ==="
  awk -v c="chr${c}" '$1==c' "${REGIONS}" > "${OUTDIR}/per_chrom/regions_chr${c}.bed"

  if [ ! -s "${OUTDIR}/per_chrom/regions_chr${c}.bed" ]; then
    echo "no variants on chr${c}, skipping"
    continue
  fi

  bcftools query \
    -f "${FMT}" \
    -R "${OUTDIR}/per_chrom/regions_chr${c}.bed" \
    "${GNOMAD_BASE}/${GNOMAD_PREFIX}.chr${c}.vcf.bgz" \
    > "${OUTDIR}/per_chrom/gnomad_chr${c}.tsv"

  wc -l "${OUTDIR}/per_chrom/gnomad_chr${c}.tsv"
done

# ---- combine -----------------------------------------------------------------
echo -e "chrom\tpos\tref\talt\tgnomAD_AF_global\tgnomAD_AC\tgnomAD_AN\tgnomAD_AF_grpmax\tgnomAD_AF_nfe\tgnomAD_AF_afr\tgnomAD_AF_eas\tgnomAD_AF_sas\tgnomAD_AF_amr" \
  > "${OUTDIR}/gnomad_af.tsv"
cat "${OUTDIR}"/per_chrom/gnomad_chr*.tsv >> "${OUTDIR}/gnomad_af.tsv"

echo "done: ${OUTDIR}/gnomad_af.tsv"
wc -l "${OUTDIR}/gnomad_af.tsv"


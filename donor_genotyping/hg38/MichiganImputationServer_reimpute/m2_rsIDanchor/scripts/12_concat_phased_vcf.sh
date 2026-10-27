#!/bin/bash
# shuo.shan@umassmed.edu
#BSUB -J concat_phased_vcfs
#BSUB -R "rusage[mem=8000]"
#BSUB -o step12_concatPhasedVCF.out
#BSUB -e step12_concatPhasedVCF.err
#BSUB -q short
#BSUB -W 2:00
#BSUB -n 8

# Concatenate the per-chromosome phased VCFs from step 11 into one genome-wide file.
# (The chromBPNet personal-genome pipeline reads the per-chromosome files directly;
#  this file is only for other uses.)

# Load necessary modules
module load bcftools
module load htslib

DIR=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor
VCF_DIR=${DIR}/MIS_results/combined_output/final_output/phased
OUT_DIR=${DIR}/MIS_results/combined_output/final_output
cd ${OUT_DIR}

# ----------- concatenate phased genotype data
# Chromosome list (space-separated). Step 11 phases autosomes only.
CHROM_LIST="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"

# VCF filename pattern (must include "{CHR}" placeholder).
PHASED_PATTERN="chr{CHR}_imputed_combined_hg38_rsID_filtered_phased.vcf.gz"

# Output file base
PHASED_OUT_BASE=skineQTL_imputed_hg38_filtered_phased

# ---- Function: Expand file list for bcftools concat ----
build_vcf_list() {
  local pattern="$1"
  local list=""
  for CHR in $CHROM_LIST; do
          file="${VCF_DIR}/$(echo ${pattern} | sed "s/{CHR}/${CHR}/")"
          if [[ -f ${file} ]]; then
                  list+="$file "
          else
                  echo "Warning: Missing file for chr$CHR: $file" >&2
          fi
  done
  echo "$list"
}

# ---- Concatenate phased VCFs ----
echo "Concatenating phased VCFs..."; date
PHASED_VCF_LIST=$(build_vcf_list ${PHASED_PATTERN})
PHASED_OUT=${OUT_DIR}/${PHASED_OUT_BASE}.vcf.gz

# chromosomes are in order and do not overlap -> plain concat, written straight to bgzip
bcftools concat --threads 8 -Oz -o ${PHASED_OUT} ${PHASED_VCF_LIST}
tabix -f ${PHASED_OUT}

n=$(bcftools view -H ${PHASED_OUT} | wc -l | cut -d' ' -f1)
echo "Done: ${PHASED_OUT} (${n} sites)"; date


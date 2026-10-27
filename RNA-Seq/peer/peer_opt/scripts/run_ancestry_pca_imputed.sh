#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=64000]
#BSUB -q long
#BSUB -W 24:00
#BSUB -J pca_imputed
#BSUB -e "./pca_imputed_%J.err"
#BSUB -o "./pca_imputed_%J.out"
# =============================================================================
# Genotype ancestry / covariate PCA on the IMPUTED, filtered VCF (the same
# genotypes used in QTL modeling). Study-samples-only PCA -> covariates.
#
# The imputed VCF is already filtered (MAF>=0.08 & F_MISSING<=0.1), so no MAF/
# missingness step is needed here. It is NOT LD-pruned, though, and high-LD
# regions (MHC, inversions) would otherwise dominate the PCs -- so we LD-prune
# and drop high-LD regions before --pca.
#
# Output for QTL: imputed_pca.eigenvec (PC coords) + imputed_pca.eigenval
#   -> feed .eigenval to the scree; reformat .eigenvec with
#      format_eigenvec_for_qtl.R into the FID/donor/donor_num/PC1..PC10 layout
#      the QTL pipeline + optimizer expect.
# =============================================================================
set -euo pipefail

module load plink2 || module load plink/2.00 || true   # need plink2 for the 6M-variant VCF

dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/ancestry_PCA_rerun2026
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz

geno=${dir}/genodata
ref=${dir}/reference
out=${dir}/output
mkdir -p "$geno" "$ref" "$out" "${dir}/plink_log"
cd "$dir"

# high-LD regions to exclude (hg38); reuse if already downloaded
highld=${ref}/high-LD-regions-hg38-GRCh38.txt
if [[ ! -s "$highld" ]]; then
  wget -q https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/high-LD-regions-hg38-GRCh38.txt -O "$highld" \
    || echo "WARN: could not fetch high-LD regions; step will error if missing"
fi

## ---- 1) VCF -> plink bed: autosomes, biallelic ACGT SNPs, unique IDs ---- ##
echo "[1/4] VCF -> bed"; date
plink2 --vcf "$vcf" \
       --autosome \
       --snps-only just-acgt \
       --max-alleles 2 \
       --set-missing-var-ids '@:#' \
       --rm-dup force-first \
       --make-bed --out "${geno}/imp"
       # --maf 0.05 --geno 0.02   # already applied in the imputed VCF (MAF>=0.08, F_MISSING<=0.1)

## ---- 2) LD-prune (drop high-LD regions), r^2>0.2 in 50-variant windows ---- ##
echo "[2/4] LD prune"; date
plink2 --bfile "${geno}/imp" \
       --exclude range "$highld" \
       --indep-pairwise 50 5 0.2 \
       --out "${geno}/imp_prune"

## ---- 3) PCA on the pruned common-variant set (study samples only) ---- ##
echo "[3/4] PCA"; date
plink2 --bfile "${geno}/imp" \
       --extract "${geno}/imp_prune.prune.in" \
       --pca 10 \
       --out "${out}/imputed_pca"
# Produces:
#   ${out}/imputed_pca.eigenvec  -> per-sample PC1..PC10 (unit-norm eigenvectors)
#   ${out}/imputed_pca.eigenval  -> the 10 eigenvalues == variance-explained for
#                                   the SCREE PLOT (send this file for the scree)

## ---- 4) Relatedness check (settle F95/F93 etc.) ---- ##
echo "[4/4] relatedness (KING)"; date
plink2 --bfile "${geno}/imp" \
       --extract "${geno}/imp_prune.prune.in" \
       --make-king-table \
       --king-table-filter 0.0884 \
       --out "${out}/relatedness"
# ${out}/relatedness.kin0 lists pairs with kinship > 0.0884 (2nd-degree or closer;
# ~0.354 = duplicate/MZ, ~0.25 = 1st-degree). To auto-drop one of each related
# pair before a final PCA:  plink2 ... --king-cutoff 0.0884 --out ${out}/unrelated
echo "done"; date

echo "Next: reformat for QTL ->"
echo "  Rscript format_eigenvec_for_qtl.R --eigenvec ${out}/imputed_pca.eigenvec --out ${out}/genotype_PCs_for_modeling_imputed.txt"

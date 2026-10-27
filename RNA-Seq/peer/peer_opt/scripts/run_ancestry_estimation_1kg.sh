#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=96000]
#BSUB -q long
#BSUB -W 24:00
#BSUB -J ancestry1kg
#BSUB -e "./ancestry1kg_%J.err"
#BSUB -o "./ancestry1kg_%J.out"
# =============================================================================
# Ancestry estimation for the study donors by projecting the IMPUTED genotypes
# onto the 1000 Genomes (GRCh38) reference (HapMap replacement).
#
#   study (imputed hg38)  +  1000G GRCh38 reference  -> combined PCA
#   -> assign_ancestry.R uses the 1000G SuperPop labels to label each donor.
#
# Both datasets are GRCh38, so no build liftover is needed -- we just normalize
# chromosome names, key variants on chr:pos, drop strand-ambiguous (AT/CG) SNPs,
# flip the rest, and merge. Study-vs-reference labeling only; the COVARIATE PCA
# for QTL is the separate study-only run_ancestry_pca_imputed.sh.
# =============================================================================
set -euo pipefail
module load plink2 || true
module load plink/1.90b6.27 || module load plink || true   # bmerge needs plink1.9
module load bcftools || true

dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/ancestry_PCA_rerun2026
studyvcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz
SCRIPTDIR=${dir}   # where assign_ancestry.R lives (edit if you keep scripts elsewhere)

geno=${dir}/genodata
refdir=${dir}/reference_1kg
out=${dir}/ancestry_output
mkdir -p "$geno" "$refdir" "$out" "${dir}/plink_log"
cd "$dir"

## ---- 0) Obtain the 1000 Genomes GRCh38 reference (plink2 pgen) ------------ ##
# Downloaded once. VERIFY current URLs at https://www.cog-genomics.org/plink/2.0/resources
# (the "1000 Genomes phase 3" GRCh38 pgen/pvar/psam). If you already have a
# reference fileset, drop ref.bed/bim/fam into $refdir and skip this block.
if [[ ! -s "${refdir}/ref.bed" ]]; then
  echo "[0] fetching 1000G GRCh38 reference"; date
  cd "$refdir"
  wget -q "https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1"  -O all_hg38.pgen.zst
  wget -q "https://www.dropbox.com/s/ngbo2xm5ojw9koy/all_hg38.pvar.zst?dl=1"  -O all_hg38.pvar.zst
  wget -q "https://www.dropbox.com/s/2e87z6nc4qexjjm/hg38_corrected.psam?dl=1" -O all_hg38.psam
  plink2 --zst-decompress all_hg38.pgen.zst > all_hg38.pgen
  # 1000G reference: autosomal biallelic ACGT SNPs, common (MAF>=0.05), chr-prefixed,
  # variant IDs = chr:pos so they align with the study by position.
  plink2 --pgen all_hg38.pgen --pvar all_hg38.pvar.zst --psam all_hg38.psam \
         --allow-extra-chr --autosome --snps-only just-acgt --max-alleles 2 \
         --maf 0.05 --output-chr chrM --set-all-var-ids '@:#' --rm-dup force-first \
         --make-bed --out "${refdir}/ref"
  cd "$dir"
fi
# reference superpopulation labels: IID + SuperPop  (from the 1000G psam)
awk 'NR==1{for(i=1;i<=NF;i++){if($i=="#IID"||$i=="IID")id=i; if($i=="SuperPop")sp=i} next}
     {print $id"\t"$sp}' "${refdir}/all_hg38.psam" > "${out}/ref_superpop.txt" 2>/dev/null || \
awk 'NR>1{print $1"\t"$5}' "${refdir}/all_hg38.psam" > "${out}/ref_superpop.txt"

## ---- 1) Study (imputed) -> bed, chr:pos IDs, drop AT/CG ambiguous --------- ##
echo "[1] study -> bed"; date
plink2 --vcf "$studyvcf" --autosome --snps-only just-acgt --max-alleles 2 \
       --output-chr chrM --set-all-var-ids '@:#' --rm-dup force-first \
       --make-bed --out "${geno}/study"
# strand-ambiguous SNPs must go before a reference merge (can't resolve flips)
awk 'BEGIN{OFS="\t"} ($5$6=="GC"||$5$6=="CG"||$5$6=="AT"||$5$6=="TA"){print $2}' \
    "${geno}/study.bim" > "${geno}/study.ac_gt_snps"
plink2 --bfile "${geno}/study" --exclude "${geno}/study.ac_gt_snps" --make-bed --out "${geno}/study.noambig"

## ---- 2) Common variants + LD-prune on the study ------------------------- ##
echo "[2] shared variants + prune"; date
awk '{print $2}' "${refdir}/ref.bim"        | sort > "${geno}/ref.ids"
awk '{print $2}' "${geno}/study.noambig.bim"| sort > "${geno}/study.ids"
comm -12 "${geno}/ref.ids" "${geno}/study.ids" > "${geno}/shared.ids"
plink2 --bfile "${geno}/study.noambig" --extract "${geno}/shared.ids" \
       --indep-pairwise 50 5 0.2 --out "${geno}/study.prune"

plink2 --bfile "${geno}/study.noambig" --extract "${geno}/study.prune.prune.in" --make-bed --out "${geno}/study.pruned"
plink2 --bfile "${refdir}/ref"         --extract "${geno}/study.prune.prune.in" --make-bed --out "${geno}/ref.pruned"

## ---- 3) Resolve strand flips, then merge (plink1.9) --------------------- ##
echo "[3] merge study + reference"; date
# first attempt -> any strand problems land in .missnp
plink --bfile "${geno}/study.pruned" --bmerge "${geno}/ref.pruned" \
      --make-bed --out "${geno}/merge_try1" || true
if [[ -s "${geno}/merge_try1-merge.missnp" ]]; then
  plink --bfile "${geno}/ref.pruned" --flip "${geno}/merge_try1-merge.missnp" \
        --make-bed --out "${geno}/ref.flipped"
  plink --bfile "${geno}/study.pruned" --bmerge "${geno}/ref.flipped" \
        --make-bed --out "${geno}/merge_try2" || true
  # anything still mismatched after flipping = true position/allele conflict: drop it
  if [[ -s "${geno}/merge_try2-merge.missnp" ]]; then
    plink --bfile "${geno}/ref.flipped" --exclude "${geno}/merge_try2-merge.missnp" --make-bed --out "${geno}/ref.final"
    plink --bfile "${geno}/study.pruned" --exclude "${geno}/merge_try2-merge.missnp" --make-bed --out "${geno}/study.final"
  else
    cp "${geno}/ref.flipped.bed" "${geno}/ref.final.bed"; cp "${geno}/ref.flipped.bim" "${geno}/ref.final.bim"; cp "${geno}/ref.flipped.fam" "${geno}/ref.final.fam"
    cp "${geno}/study.pruned.bed" "${geno}/study.final.bed"; cp "${geno}/study.pruned.bim" "${geno}/study.final.bim"; cp "${geno}/study.pruned.fam" "${geno}/study.final.fam"
  fi
else
  cp "${geno}/ref.pruned.bed" "${geno}/ref.final.bed"; cp "${geno}/ref.pruned.bim" "${geno}/ref.final.bim"; cp "${geno}/ref.pruned.fam" "${geno}/ref.final.fam"
  cp "${geno}/study.pruned.bed" "${geno}/study.final.bed"; cp "${geno}/study.pruned.bim" "${geno}/study.final.bim"; cp "${geno}/study.pruned.fam" "${geno}/study.final.fam"
fi
plink --bfile "${geno}/study.final" --bmerge "${geno}/ref.final" --make-bed --out "${out}/combined"

## ---- 4) Combined PCA ---------------------------------------------------- ##
echo "[4] combined PCA"; date
plink2 --bfile "${out}/combined" --pca 10 --out "${out}/combined_pca"
# study sample IDs (to split ref vs study downstream)
awk '{print $2}' "${geno}/study.final.fam" > "${out}/study_iids.txt"

## ---- 5) Assign ancestry (k-NN on reference SuperPop) -------------------- ##
echo "[5] ancestry assignment"; date
Rscript "${SCRIPTDIR}/assign_ancestry.R" \
  --eigenvec "${out}/combined_pca.eigenvec" \
  --ref-labels "${out}/ref_superpop.txt" \
  --study-iids "${out}/study_iids.txt" \
  --out "${out}/ancestry_estimation.txt" \
  --plot "${out}/ancestry_PCA_1kg.pdf" || echo "assign_ancestry.R failed; check paths"

echo "done"; date
echo "Ancestry table: ${out}/ancestry_estimation.txt"

#!/bin/bash
# fetch genotype of all tested SNPs of a given gene

SCRIPTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/script
OUTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats
module load bcftools

# ---------- fetch input arguments ------------------- #
g=$1   #g=ITGA1
letter=${g:0:1}
mkdir -p ${OUTDIR}/${letter}/${g}
cd ${OUTDIR}/${letter}/${g}
echo "starting the script for ${g}"; date

# get SNP genotype
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz
vcf_subset=snp.vcf.gz
gene_snp_pairs=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/SNPs_near_TSS.bed
samples=samples.txt
header=header
body=body
geno=temp.snp.txt
geno_num=genotype.txt

echo "==== Subsetting genotype VCF for tested SNPs of ${g} ===="; date
# 1) collect SNP IDs for all tested SNPs of this gene
file="${OUTDIR}/${letter}/${g}/compiled_modeling_stats_${g}.txt"
if [ -f "${OUTDIR}/${letter}/compiled_modeling_stats_${g}.txt" ]; then
    mv "${OUTDIR}/${letter}/compiled_modeling_stats_${g}.txt" "$file"
fi

awk 'NR>1{print $1}' ${file} | sort -u > snplist.txt

# 2) extract VCF records
bcftools view --include ID==@snplist.txt ${vcf} -Oz -o ${vcf_subset}

# 3) clean sample names
bcftools query -l "$vcf_subset" | cut -d'_' -f1 | sed 's/F0/F/g' | sed 's/skineQTL-//g' > "$samples"

# 4) build header
printf "CHROM\tPOS\tID\tREF\tALT\t" > "$header"
paste -sd '\t' "$samples" >> "$header"

# 5) extract GT matrix
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' "$vcf_subset" > "$body"
cat "$header" "$body" > "$geno"

# 6) convert genotypes to numeric
sed -E 's/0\/0/0/g; s/0\/1/1/g; s/1\/0/1/g; s/1\/1/2/g; s/\.\/\./NA/g' "$geno" > "$geno_num"

# 7) cleanup intermediate files (leave only numeric genotype)
rm snplist.txt "$vcf_subset" "$samples" "$header" "$body" "$geno"

echo "SNP genotype written to ${OUTDIR}/${letter}/${g}/genotype.txt"

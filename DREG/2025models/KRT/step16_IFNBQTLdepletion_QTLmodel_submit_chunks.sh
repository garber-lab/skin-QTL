#!/bin/bash
#BSUB -n 5
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -W 72:00
#BSUB -q long
#BSUB -J QTL.KRT[1-135]%70
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/KRT/IFNBeGene_depletion_issue/logs/qtl_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/KRT/IFNBeGene_depletion_issue/logs/qtl_%J_%I.err"

sleep $((RANDOM % 20)) # reduces heavy I/O thundering-herd effect in parallel jobs
module load bcftools
module load htslib

celltype=KRT
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${celltype}/IFNBeGene_depletion_issue
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz

# Select chunk file safely
chunk=$(ls ${DIR}/chunks/pairs_chunk_*.tsv | sort -V | sed -n ${LSB_JOBINDEX}p)
#chunk=${DIR}/chunks/pairs_chunk_096.tsv

# Extract chunk base name and chunk index (000, 001, ...)
chunk_base=$(basename "$chunk" .tsv)
chunk_id=$(echo "$chunk_base" | sed 's/pairs_chunk_//')

# Define output file path
out="${DIR}/results/result_${chunk_id}.tsv"

# Define intermediate filenames (never overwrite $chunk)
snplist=${DIR}/chunks/temp.${chunk_base}.snps.txt
vcf_subset=${DIR}/chunks/temp.${chunk_base}.vcf.gz
samples=${DIR}/chunks/temp.${chunk_base}.samples.txt
header=${DIR}/chunks/temp.${chunk_base}.header.txt
body=${DIR}/chunks/temp.${chunk_base}.body.tsv
geno=${DIR}/chunks/temp.${chunk_base}.genotype.tsv
geno_num=${DIR}/chunks/genotype_${chunk_base}.tsv


# -------- Subset genotype VCF file to 100K --------------------- #
echo "Subsetting genotype VCF file for ${chunk}";date
cd ${DIR}/chunks

echo "==== Subsetting genotype VCF for $chunk ===="; date

# 1) collect SNP IDs from this chunk
awk 'NR>1 {print $11}' "$chunk" | sort -u > "$snplist"

# 2) extract VCF records
bcftools view --threads 5 --include ID==@"$snplist" "$vcf" -Oz -o "$vcf_subset"

# 3) clean sample names
bcftools query -l ${vcf_subset} | cut -d'_' -f1 | sed 's/skineQTL-//' | sed 's/F0/F/' > ${samples}

# 4) build header
printf "CHROM\tPOS\tID\tREF\tALT\t" > "$header"
paste -sd '\t' "$samples" >> "$header"

# 5) extract GT matrix
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' "$vcf_subset" > "$body"
cat "$header" "$body" > "$geno"

# 6) convert genotypes to numeric
sed -E 's/0\/0/0/g; s/0\/1/1/g; s/1\/0/1/g; s/1\/1/2/g; s/\.\/\./NA/g' "$geno" > "$geno_num"

# 7) cleanup intermediate files (leave only numeric genotype)
rm "$snplist" "$vcf_subset" "$samples" "$header" "$body" "$geno"


# -------- Set-up Singularity and Dependencies for Rscript ------ #
echo "Fitting QTL model for ${chunk}";date
cd ${DIR}

# Tell R to use my library inside the container:
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${celltype}/step16_IFNBQTLdepletion_QTLmodel_run_QTL_chunk.R ${celltype} 10 2 ${chunk_id}

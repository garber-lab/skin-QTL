#!/bin/bash
# shuo.shan@umassmed.edu
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -W 08:00
#BSUB -q long
#BSUB -J QTL
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step0_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step0_%J_%I.err"

# ----------------------------------------- #
# -------- Global Variables        -------- #
# ----------------------------------------- #
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models
SCRIPTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/scripts
metagene=/pi/manuel.garber-umw/human/skin/eQTLs/literature/metaidname.txt
genemodel=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/data/gencode_v34_gene_TSS.bed
genelist=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/data/expressed_genes_all.txt
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz
pairs=${DIR}/data/SNPs_near_TSS.bed
CHUNK_DIR=${DIR}/data/chunks

cd ${DIR}
mkdir -p ${CHUNK_DIR}
mkdir -p ${DIR}/data
mkdir -p ${DIR}/logs
mkdir -p ${DIR}/results

# ----------------------------------------- #
# -------- Fetch expressed gene list ------ #
# ----------------------------------------- #
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
for ct in MEL KRT FRB;do
	for cond in PBS IFNG IFNB TNF;do
		singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
			Rscript ${SCRIPTDIR}/step0.1_pick_expressed_genes.R ${ct} ${cond} expressed
	done

	for cond in IFNG IFNB TNF;do
		singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
                        Rscript ${SCRIPTDIR}/step0.1_pick_expressed_genes.R ${ct} ${cond} induced
	done
done

cd ${DIR}/data
for f in gene_list_*.txt;do
	cat ${f} >> ${DIR}/data/temp.gene_list.txt
done
sort -u ${DIR}/data/temp.gene_list.txt > ${genelist}
rm ${DIR}/data/temp.gene_list.txt

# ----------------------------------------- #
# -------- Build gene TSS ----------------- #
# ----------------------------------------- #

awk 'BEGIN{FS="\t"; OFS="\t"}
$3=="gene" {
  n = split($9, a, ";")
  gene_id=""; gene_name=""
  for (i=1; i<=n; i++) {
    if (a[i] ~ /gene_id/)   { s=a[i]; gsub(/^[ ]*gene_id "/,"",s); gsub(/"$/,"",s); gene_id=s }
    if (a[i] ~ /gene_name/) { s=a[i]; gsub(/^[ ]*gene_name "/,"",s); gsub(/"$/,"",s); gene_name=s }
  }
  chr=$1; start=$4; end=$5; strand=$7
  if (strand=="+") { tss0=start-1; tss1=start } else { tss0=end-1; tss1=end }
  print chr, tss0, tss1, gene_name, 0, strand, gene_id
}' /pi/manuel.garber-umw/home/yuqing.wang-umw/inDrop/Skin/diablo/build/hg38_v34/gencode.v34.primary_assembly.annotation.gtf | \
        sort -k1,1 -k2,2n > ${genemodel}

# ----------------------------------------- #
# -------- Make SNP-gene pair list -------- #
# ----------------------------------------- #
module load bcftools
module load bedtools
module load htslib

# Convert VCF → BED format (1-based → 0-based)
bcftools view -H ${vcf} | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' > ${DIR}/data/all_genotyped_SNPs.bed

# Print header
echo -e "gene_chr\tgene_start\tgene_end\tgene_name\t0\tstrand\tgene_id\tSNP_chr\tSNP_start\tSNP_end\tSNP_ID" > ${pairs}

# Get every SNP within +/- 500kb of a gene's TSS
bedtools window -w 500000 -a ${genemodel} -b ${DIR}/data/all_genotyped_SNPs.bed > ${pairs}.temp

# Filter to the list of genes (expressed in any context)
awk 'NR==FNR{genes[$1]; next} $4 in genes' ${genelist} ${pairs}.temp >> ${pairs}

# clean-up
rm ${pairs}.temp

# ----------------------------------------- #
# ------ chunk genes -------------- ####
# ----------------------------------------- #
# chunk gene list so that each chunk correspond to <= 50K gene_snp pairs
echo "Chunking gene list with gene:snp cap of 50K per chunk..."

# 1) make a stable unique gene list
all_genes=${CHUNK_DIR}/all_genes.txt
tail -n +2 ${pairs} \
  | awk -F'\t' '{print $4}' \
  | awk '!seen[$0]++' \
  | shuf \
  > ${all_genes}

# 2) split gene list into max of 100,000 variant-gene pair chunks
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step0.2_compile_gene_chunks.R 
echo "Num gene chunks: $(ls -1 ${CHUNK_DIR}/gene_chunk_* | wc -l)"
# 14599 genes into 314 chunks | genes/chunk: min=1 median=47 max=50

# 3) for each gene chunk, extract all matching rows from the big pairs file (SNPs are bed 0-based coordinated now)
#    (add header to each output)
header=$(head -n1 ${pairs})

for gfile in ${CHUNK_DIR}/gene_chunk_*; do
  chunk_id=$(basename ${gfile} | sed 's/gene_chunk_//')
  out_pairs=${CHUNK_DIR}/pairs_chunk_${chunk_id}.tsv

  echo "Building ${out_pairs} from genes in ${gfile}"; date

  awk -F'\t' -v OFS='\t' -v header="${header}" '
    NR==FNR { genes[$1]=1; next }     # read gene list into hash
    FNR==1 { print header; next }      # at first line of pairs file: print header, skip it
    ($4 in genes) { print }            # print matching pair rows
  ' "${gfile}" "${pairs}" > "${out_pairs}"
done

echo "Done. Gene-chunked pair files in ${CHUNK_DIR}/"


# ----------------------------------------- #
# ------ chunk VST per gene chunk --------- ####
# ----------------------------------------- #
# Write a per-chunk VST subset (vst_chunk_<id>.tsv) matching each gene_chunk_<id>,
# so step1 reads a tiny file instead of the full ~409MB VST (much less memory per
# job, no repeated big read). Does NOT touch the gene/genotype chunks. Idempotent.
VST_FILE=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt
echo "Chunking VST to match gene chunks..."; date
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step0.4_chunk_vst.R ${VST_FILE} ${CHUNK_DIR}


# -------------------------------------------------- #
# ----- Compile Genotype table for each chunk ------ #
# -------------------------------------------------- #
module load bcftools
module load htslib

# Select chunk file safely
for chunk in ${CHUNK_DIR}/pairs_chunk_*.tsv; do

	# Extract chunk base name and chunk index (000, 001, ...)
	chunk_id=$(basename ${chunk} | sed 's/pairs_chunk_//' | sed 's/\.tsv//g' )
	chunk_base=pairs_chunk_${chunk_id}
	
	# Define intermediate filenames (never overwrite $chunk)
	snplist=${CHUNK_DIR}/temp.${chunk_base}.snps.txt
	vcf_subset=${CHUNK_DIR}/temp.${chunk_base}.vcf.gz
	samples=${CHUNK_DIR}/temp.${chunk_base}.samples.txt
	header=${CHUNK_DIR}/temp.${chunk_base}.header.txt
	body=${CHUNK_DIR}/temp.${chunk_base}.body.tsv
	geno=${CHUNK_DIR}/temp.${chunk_base}.genotype.tsv
	geno_num=${CHUNK_DIR}/genotype_${chunk_base}.tsv
	
	
	# -------- Subset genotype VCF file to 100K --------------------- #
	echo "Subsetting genotype VCF file for ${chunk}";date
	cd ${CHUNK_DIR}
	
	
	# 1) collect SNP IDs from this chunk
	awk 'NR>1 {print $11}' "$chunk" | sort -u > "$snplist"
	
	# 2) extract VCF records
	bcftools view --threads 8 --include ID==@"$snplist" "$vcf" -Oz -o "$vcf_subset"
	
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

done




# generate the right job-array
echo "#BSUB -J QTL[1-314]%70"
# #BSUB -J QTL[1-314]%70 # create 314 jobs but run only 70 at a time concurrently, queue the rest

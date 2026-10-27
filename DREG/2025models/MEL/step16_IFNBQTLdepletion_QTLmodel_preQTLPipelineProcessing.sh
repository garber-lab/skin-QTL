#!/bin/bash
# shuo.shan@umassmed.edu
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -W 08:00
#BSUB -q long
#BSUB -J QTL.MEL
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/step0_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/step0_%J_%I.err"

# ----------------------------------------- #
# -------- Global Variables        -------- #
# ----------------------------------------- #
ct=MEL
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}/IFNBeGene_depletion_issue
metagene=/pi/manuel.garber-umw/human/skin/eQTLs/literature/metaidname.txt
genemodel=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/gencode_v34_gene_TSS.bed
genelist=${DIR}/data/expressed_and_induced_genes_${ct}_IFNB_or_IFNG.txt
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz
pairs=${DIR}/data/SNPs_near_TSS_induced.bed
CHUNK_DIR=${DIR}/chunks

cd ${DIR}

# ----------------------------------------- #
# -------- Build gene TSS ----------------- #

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

# Filter to the list of genes (expressed and induced in IFNB or IFNG)
awk 'NR==FNR{genes[$1]; next} $4 in genes' ${genelist} ${pairs}.temp >> ${pairs}

# clean-up
rm ${pairs}.temp

# ------ chunk genes -------------- ####
# chunk gene list so that each chunk correspond to <= 50K gene_snp pairs
echo "Chunking gene list with gene:snp cap of 50K per chunk..."

# 1) make a stable unique gene list
all_genes=${CHUNK_DIR}/all_genes.txt
tail -n +2 ${pairs} \
  | awk -F'\t' '{print $4}' \
  | awk '!seen[$0]++' \
  | shuf \
  > ${all_genes}

# 2) split gene list into max of 50,000 variant-gene pair chunks, min 11 gene, mediam 20 genes, max 20 genes
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
        Rscript /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}/step16_IFNBQTLdepletion_QTLmodel_decide_gene_chunks.R ${ct}
echo "Num gene chunks: $(ls -1 ${CHUNK_DIR}/gene_chunk_* | wc -l)"

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

# -------------------

# detect total number of chunks automatically
NCHUNKS=$(ls ${DIR}/chunks/pairs_chunk_*.tsv | wc -l)
echo "Total chunks: ${NCHUNKS}"
# 83 chunks

# generate the right job-array
echo "#BSUB -J QTL[1-${NCHUNKS}]%70"
echo "#BSUB -J QTL[1-83]%70"
# #BSUB -J QTL[1-273]%70 # create 273 jobs but run only 70 at a time concurrently, queue the rest

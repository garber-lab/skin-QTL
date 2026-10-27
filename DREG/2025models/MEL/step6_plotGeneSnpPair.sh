#!/bin/bash
# make a PDF for gene snp pair of interest given the SNP is in the gene's cis-window.

ct=MEL
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
#mkdir -p ${DIR}/plots

# ---------- fetch input arguments ------------------- #
g=$1   #g=ITGA1
snp=$2 #snp=rs2548496
letter=${g:0:1}
mkdir -p ${DIR}/plots/${letter}/${g}
mkdir -p ${DIR}/plots/${letter}/${g}/${snp}
cd ${DIR}/plots/${letter}/${g}/${snp}
echo "starting the script for ${g}:${snp}"; date

## ---------- create gene chunk lookup table ---------- #
#cd ${DIR}/data/chunk
#echo -e "chunk\tgene" > ${DIR}/data/gene_chunk_dict.txt
#for f in gene_chunk_*;do
#	id=$(echo ${f} | cut -d"_" -f3)
#	echo ${f} " ID is: " ${id}
#	awk -v id=${id} '{OFS="\t"}{print id,$0}' ${f} >> ${DIR}/data/gene_chunk_dict.txt
#done
#
## ---------- create gene chunk lookup table from old run ----------- #
#cd ${DIR}/results
#echo -e "chunk\tgene" > ${DIR}/data/gene_chunk_dict_old.txt
#for f in result_*.tsv;do
#	id=$(echo ${f} | cut -d"_" -f2 | cut -d"." -f1)
#	echo ${f} " ID is: " ${id}
#	awk -F "\t" 'NR>1{print $2}' result_${id}.tsv | sort -u > temp_gene_${id}.txt
#	awk -v id=${id} '{OFS=FS="\t"}{print id,$1}' temp_gene_${id}.txt >> ${DIR}/data/gene_chunk_dict_old.txt
#done

# ---------- get gene, SNP, metadata, and modeling stats info ----------- #
# for quick lookup, fetch chunk ID of gene
chunk_dict=${DIR}/data/gene_chunk_dict.txt
chunk_id=$(awk -v gene=${g} '{if ($2==gene) print $1}' ${chunk_dict})

# get gene CPM
cpm_all_f=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt
head -1 ${cpm_all_f} > ${DIR}/plots/${letter}/${g}/${snp}/cpm.txt
awk -v gene=${g} '{OFS=FS="\t"}{if ($NF==gene) print $0}' ${cpm_all_f} >> ${DIR}/plots/${letter}/${g}/${snp}/cpm.txt
echo "got gene CPM"

# get SNP genotype
genotype_chunk_f=${DIR}/data/chunk/genotype_pair_chunk_${chunk_id}.txt
head -1 ${genotype_chunk_f} > ${DIR}/plots/${letter}/${g}/${snp}/genotype.txt
awk -v snp=${snp} '{OFS=FS="\t"}{if ($3==snp) print $0}' ${genotype_chunk_f} >> ${DIR}/plots/${letter}/${g}/${snp}/genotype.txt 
echo "got SNP genotype"

# get model stats
# get all chunk ids for this gene (may be 1 or many)
chunk_id_old=$(awk -v gene="$g" '$2==gene {print $1}' ${DIR}/data/gene_chunk_dict_old.txt)

out="${DIR}/plots/${letter}/${g}/${snp}/modelstats.txt"
mkdir -p "$(dirname "$out")"

# write header from first chunk
first_chunk=$(echo "$chunk_id_old" | head -n1)
echo -e "celltype\t$(head -1 ${DIR}/results/result_${first_chunk}.tsv)" > "$out"

# loop through all chunks, append match if exists
for cid in $chunk_id_old; do
  awk -v snp="$snp" -v gene="$g" -v ct="$ct" '
    BEGIN{FS=OFS="\t"}
    NR>1 && $1==snp && $2==gene {print ct, $0}
  ' ${DIR}/results/result_${cid}.tsv >> "$out"
done

echo "got modeling stats"

# ------------ call Rscript to plot --------------- #
echo "plotting now..."
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${DIR}/step6_make_plot_for_pair_in_window.R ${ct} ${g} ${snp}

cd ${DIR}/plots
rm -r ${DIR}/plots/${letter}/${g}/${snp}
rm -r ${DIR}/plots/${letter}/${g}

echo "done!"; date











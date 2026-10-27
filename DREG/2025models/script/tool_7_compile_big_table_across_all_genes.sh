#!/bin/bash
# workflow to run FIMO on all shortlisted chromatin-disrupting variants for all genes (output from tool8)

SCRIPTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/script
OUTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats

# ---- fetch all variants and their genes ----
gene_file=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_genes.bed
gene_list=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_genes.txt
cut -f4 ${gene_file} | awk 'NR>1' | sort -u > ${gene_list}

chromBPNet_variant_score_QTL_masterfile=${OUTDIR}/chromBPNet_variant_score_QTL_masterfile.txt
rm -f ${chromBPNet_variant_score_QTL_masterfile}

head -1 /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats/G/GBP3/big_table.txt > ${chromBPNet_variant_score_QTL_masterfile}

while IFS= read -r g; do

    letter=${g:0:1}
    big_table_file=${OUTDIR}/${letter}/${g}/big_table.txt
    awk 'NR>1{print $0}' ${big_table_file} >> ${chromBPNet_variant_score_QTL_masterfile}

    echo "[DONE] ${g}"

done < ${gene_list}

echo "Master file written to: ${chromBPNet_variant_score_QTL_masterfile}"


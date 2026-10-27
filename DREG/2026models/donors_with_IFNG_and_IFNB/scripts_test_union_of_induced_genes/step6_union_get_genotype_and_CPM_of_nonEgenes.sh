#!/bin/bash
# UNION variant of step6_get_genotype_and_CPM_of_nonEgenes.sh: extract compact
# genotype + CPM tables for the matched NON-eGene control set (UNION experiment).
# Inputs/outputs live in <ROOT>/results_test_union_of_induced_genes.
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4000]"
#BSUB -W 02:00
#BSUB -q short
#BSUB -J getGenoCPM_ctrl_union
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/step6_union_getGenoCPM_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/step6_union_getGenoCPM_%J.err"

set -e
ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB
CHUNK_DIR=${ROOT}/data/chunks
RES=${ROOT}/results_test_union_of_induced_genes
CPM_ALL=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt

OUT_GENO=${RES}/genotype_of_control_nonEgene_tested_variants.txt
OUT_CPM=${RES}/CPM_of_control_nonEgenes.txt

# compiled step6 control columns: celltype(1) condition(2) QTLtype(3) snp(4) gene(5) ...
IN_GLOB="${RES}/QTL_output_of_all_tested_variants_of_control_nonEgenes_*.txt"
if ! ls ${IN_GLOB} >/dev/null 2>&1; then
  echo "ERROR: no union step6 compiled control files matching ${IN_GLOB} -- run step6_union_match first." >&2
  exit 1
fi

tmp=$(mktemp -d)
awk -F'\t' 'FNR>1{print $5}' ${IN_GLOB} | sort -u > ${tmp}/genes.txt
awk -F'\t' 'FNR>1{print $4}' ${IN_GLOB} | sort -u > ${tmp}/snps.txt
echo "control genes: $(wc -l < ${tmp}/genes.txt) | tested variants: $(wc -l < ${tmp}/snps.txt)"; date

# --- CPM rows for control genes (final_gene is the LAST column) ---
head -1 ${CPM_ALL} > ${OUT_CPM}
awk -F'\t' 'NR==FNR{g[$1];next} ($NF in g)' ${tmp}/genes.txt <(tail -n +2 ${CPM_ALL}) \
  | awk -F'\t' '!seen[$NF]++' >> ${OUT_CPM}
echo "wrote CPM: $(($(wc -l < ${OUT_CPM}) - 1)) genes -> ${OUT_CPM}"

# --- genotype rows for tested variants (scan chunk genotype tables, dedup by ID) ---
first=$(ls ${CHUNK_DIR}/genotype_pairs_chunk_*.tsv | head -1)
head -1 ${first} > ${OUT_GENO}
for f in ${CHUNK_DIR}/genotype_pairs_chunk_*.tsv; do
  awk -F'\t' 'NR==FNR{s[$1];next} FNR>1 && ($3 in s)' ${tmp}/snps.txt ${f}
done | awk -F'\t' '!seen[$3]++' >> ${OUT_GENO}
echo "wrote genotype: $(($(wc -l < ${OUT_GENO}) - 1)) variants -> ${OUT_GENO}"; date

rm -rf ${tmp}

#!/bin/bash
# Plot one gene:SNP pair on demand.
#   Usage: bash step7_make_plot_single_pair.sh <GENE> <SNP_ID>
# Assumes the gene and SNP live in the same chunk (SNP is in the gene's cis-window),
# so the gene's chunk id locates both its subset QTL stats and (for excluded donors)
# its genotype.
#
#   - genotype (ALL genome-wide donors) <- /pi/.../2026models/data/chunks   (so excluded
#     donors render as light dots)
#   - CPM (all samples)                 <- global CPM.sampleFiltered.metaConverted
#   - QTL stats (subtitles)             <- this analysis' <ROOT>/<ct>/results_QC (LONG)
set -e
g=$1
snp=$2
[ -z "$g" ] || [ -z "$snp" ] && { echo "Usage: bash step7_make_plot_single_pair.sh <GENE> <SNP_ID>"; exit 1; }

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB
STIMS=IFNG,IFNB
GW_CHUNKS=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/data/chunks   # genome-wide: all donors
CPM_ALL=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt
SCRIPTDIR=${ROOT}/scripts
PLOTDIR=${ROOT}/plots
mkdir -p ${PLOTDIR}/temp_output

# ---- subset gene->chunk dict (idempotent) : locate the gene's chunk ----
dict=${ROOT}/data/gene_chunk_dict.txt
if [ ! -s "${dict}" ]; then
  echo -e "chunk\tgene" > ${dict}
  for f in ${ROOT}/data/chunks/gene_chunk_*; do
    id=$(basename ${f} | sed 's/gene_chunk_//')
    awk -v id=${id} '{print id"\t"$0}' ${f} >> ${dict}
  done
fi
chunk_id=$(awk -F'\t' -v gene=${g} '$2==gene{print $1; exit}' ${dict})
echo "gene ${g} -> chunk ${chunk_id:-<not in subset>}"

wd=$(mktemp -d)

# ---- CPM (global, all samples) for the gene ----
head -1 ${CPM_ALL} > ${wd}/cpm.txt
awk -F'\t' -v g=${g} '$NF==g' ${CPM_ALL} >> ${wd}/cpm.txt

# ---- genotype (all genome-wide donors) for the SNP ----
first=$(ls ${GW_CHUNKS}/genotype_pairs_chunk_*.tsv | head -1)
head -1 ${first} > ${wd}/genotype.txt
for f in ${GW_CHUNKS}/genotype_pairs_chunk_*.tsv; do
  awk -F'\t' -v s=${snp} 'FNR>1 && $3==s' ${f}
done | awk -F'\t' '!seen[$3]++' >> ${wd}/genotype.txt

# ---- LONG QTL stats from this analysis' results_QC (only the gene's chunk) ----
stats=${wd}/modelstats.txt
echo -e "celltype\tcondition\tQTLtype\tsnp\tgene\tdist\tbeta\tp\tse\tn_geno0\tn_geno1\tn_geno2\tqc_genotypes" > ${stats}
if [ -n "${chunk_id}" ]; then
  for ct in FRB KRT MEL; do
    for f in ${ROOT}/${ct}/results_QC/modeling_stats_postQC_${ct}_*_${chunk_id}.txt; do
      [ -e "$f" ] || continue
      bn=$(basename "$f" .txt); rest=${bn#modeling_stats_postQC_${ct}_}   # <cond>_<qtl>_<chunk>
      cond=$(echo ${rest} | cut -d'_' -f1); qtl=$(echo ${rest} | cut -d'_' -f2)
      awk -F'\t' -v OFS='\t' -v ct=${ct} -v cond=${cond} -v qtl=${qtl} -v g=${g} -v s=${snp} \
        'FNR>1 && $1==s && $2==g {print ct,cond,qtl,$0}' "$f" >> ${stats}
    done
  done
fi

# ---- plot ----
OUT=${PLOTDIR}/temp_output/plot_${g}_${snp}.pdf
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
  Rscript ${SCRIPTDIR}/step7_make_plot_gene_snp_pair.R \
    ${ROOT} ${g} ${snp} ${wd}/genotype.txt ${wd}/cpm.txt ${stats} ${OUT} ${STIMS}

rm -rf ${wd}
echo "made ${OUT}"

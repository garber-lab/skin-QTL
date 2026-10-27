#!/bin/bash
# Plot one gene:SNP pair on demand, GENOME-WIDE run.
#   Usage: bash step7_make_plot_single_pair.sh <GENE> <SNP_ID>
# Assumes the SNP is in the gene's cis-window (same chunk). PDF is written into an
# A-Z folder named after the gene's first character (ERAP2 -> E/, GBP3 -> G/).
set -e
g=$1
snp=$2
[ -z "$g" ] || [ -z "$snp" ] && { echo "Usage: bash step7_make_plot_single_pair.sh <GENE> <SNP_ID>"; exit 1; }

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models
STIMS=IFNG,IFNB,TNF
CHUNKS=${ROOT}/data/chunks
CPM_ALL=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt
SCRIPTDIR=${ROOT}/scripts
PLOTDIR=${ROOT}/plots

# A-Z folder from the gene's first character (non-letters -> 0-9)
letter=$(printf '%s' "${g:0:1}" | tr '[:lower:]' '[:upper:]')
case "$letter" in [A-Z]) ;; *) letter="0-9" ;; esac
OUTDIR=${PLOTDIR}/${letter}
mkdir -p ${OUTDIR}

# ---- gene->chunk dict (idempotent) ----
dict=${ROOT}/data/gene_chunk_dict.txt
if [ ! -s "${dict}" ]; then
  echo -e "chunk\tgene" > ${dict}
  for f in ${CHUNKS}/gene_chunk_*; do
    id=$(basename ${f} | sed 's/gene_chunk_//')
    awk -v id=${id} '{print id"\t"$0}' ${f} >> ${dict}
  done
fi
chunk_id=$(awk -F'\t' -v gene=${g} '$2==gene{print $1; exit}' ${dict})
[ -n "${chunk_id}" ] || { echo "ERROR: gene ${g} not found in ${dict}"; exit 1; }
echo "gene ${g} -> chunk ${chunk_id} -> folder ${letter}"

wd=$(mktemp -d)

# ---- CPM (global, all samples) for the gene ----
head -1 ${CPM_ALL} > ${wd}/cpm.txt
awk -F'\t' -v g=${g} '$NF==g' ${CPM_ALL} >> ${wd}/cpm.txt

# ---- genotype (all donors) for the SNP, from the gene's chunk ----
gfile=${CHUNKS}/genotype_pairs_chunk_${chunk_id}.tsv
head -1 ${gfile} > ${wd}/genotype.txt
awk -F'\t' -v s=${snp} 'FNR>1 && $3==s' ${gfile} >> ${wd}/genotype.txt

# ---- LONG QTL stats from results_QC (only the gene's chunk) ----
stats=${wd}/modelstats.txt
echo -e "celltype\tcondition\tQTLtype\tsnp\tgene\tdist\tbeta\tp\tse\tn_geno0\tn_geno1\tn_geno2\tqc_genotypes" > ${stats}
for ct in FRB KRT MEL; do
  for f in ${ROOT}/${ct}/results_QC/modeling_stats_postQC_${ct}_*_${chunk_id}.txt; do
    [ -e "$f" ] || continue
    bn=$(basename "$f" .txt); rest=${bn#modeling_stats_postQC_${ct}_}   # <cond>_<qtl>_<chunk>
    cond=$(echo ${rest} | cut -d'_' -f1); qtl=$(echo ${rest} | cut -d'_' -f2)
    awk -F'\t' -v OFS='\t' -v ct=${ct} -v cond=${cond} -v qtl=${qtl} -v g=${g} -v s=${snp} \
      'FNR>1 && $1==s && $2==g {print ct,cond,qtl,$0}' "$f" >> ${stats}
  done
done

# ---- plot ----
OUT=${OUTDIR}/plot_${g}_${snp}.pdf
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
  Rscript ${SCRIPTDIR}/step7_make_plot_gene_snp_pair.R \
    ${ROOT} ${g} ${snp} ${wd}/genotype.txt ${wd}/cpm.txt ${stats} ${OUT} ${STIMS}

rm -rf ${wd}
echo "made ${OUT}"

#!/bin/bash
# fetch CPM and VST expression data of a given gene

SCRIPTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/script
OUTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats
module load bedtools

# ---------- fetch input arguments ------------------- #
g=$1   #g=ITGA1
letter=${g:0:1}
mkdir -p ${OUTDIR}/${letter}/${g}
cd ${OUTDIR}/${letter}/${g}
echo "starting the script for ${g}"; date

# --- fetch genotype ---
VST_FILE="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
head -1 ${VST_FILE} > vst.txt
awk -v gene=${g} '{OFS=FS="\t"}{if ($1==gene) print $0}' ${VST_FILE} >> vst.txt
echo "VST written to ${OUTDIR}/${letter}/${g}/vst.txt"

# --- fetch CPM ---
CPM_FILE="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
    Rscript - ${CPM_FILE} ${g} ${OUTDIR}/${letter}/${g}/cpm.txt <<'RS'
args <- commandArgs(trailingOnly = TRUE)
cpm_file <- args[1]
target_gene <- args[2]
out_file <- args[3]

df <- read.delim( cpm_file, header = T, sep = "\t")
df <- df[df[[ncol(df)]] == target_gene, , drop = FALSE]
df <- df[, !(names(df) %in% c("name", "gene")), drop = FALSE]
names(df)[names(df) == "final_gene"] <- "gene"
df <- df[, c("gene", setdiff(names(df), "gene")), drop = FALSE]

write.table(df, file = out_file, sep = "\t", quote = F, row.names = F, col.names = T)
RS
echo "CPM written to ${OUTDIR}/${letter}/${g}/cpm.txt"


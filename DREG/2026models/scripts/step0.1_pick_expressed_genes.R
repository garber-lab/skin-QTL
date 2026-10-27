#!/bin/Rscript
# output gene list that's expressed or induced (expressed and induced) in the given context 
# expressed = CPM ≥ 1 in at least 20% of samples, and CPM ≥ 10 in at least 5 samples.
# induced = log2FC > 1, padj < 0.05
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

# input arguemnts
args        <- commandArgs(trailingOnly = TRUE)
celltype   <- args[1]
condition <- args[2]
type <- args[3] # "expressed or induced"

dir <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/data"
outfile= paste0(dir,"/gene_list_",celltype,"_",condition,"_",type,".txt")
# celltype="MEL";condition="IFNG"; type="expressed"

# select expressed genes for given celltype x condition of interest
path.Rscript.fetch_expressed_genes="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_expressed_genes.R"
source(path.Rscript.fetch_expressed_genes)
path.Rscript.fetch_induced_genes="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_induced_genes_log2FC1_padj0.05.R"
source(path.Rscript.fetch_induced_genes)

if (type=="expressed") {
  gene_list <- unique(fetch_expressed_genes(celltype, condition))
} else if (type == "induced"){
  ind <- unique(fetch_induced_genes(celltype, condition))
  exp <- unique(fetch_expressed_genes(celltype, condition))
  gene_list <- intersect(exp, ind)
}

write(gene_list, file = outfile)
message(paste(length(gene_list), outfile))

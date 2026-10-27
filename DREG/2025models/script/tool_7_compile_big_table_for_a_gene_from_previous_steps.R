#!/usr/bin/env Rscript
# For a gene of interest, annotate:
# 1. gene's DESeq result in all celltype x conditions
# 2. gene's expression status in all celltype x conditions
# 3. all tested SNP modeling stats
# 4. all tested SNP ATACseq overlapping status
# 5. all tested SNP chromBPNet variant scores
# 6. all tested SNP phyloP score

Sys.setenv(TZ = "America/New_York")
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(magrittr)
  library(dplyr)
})

# ---------------------- Argument Parsing ----------------------
args <- commandArgs(trailingOnly = TRUE)
this_gene <- args[1] 
#this_gene <- "AASS"
letter <- substr(this_gene, 1, 1)
workingdir = paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats/", letter,"/",this_gene)

## ---- 1. COMPILE BIG TABLE -----
# load modeling stats
modelstats <- fread(paste0(workingdir,"/compiled_modeling_stats_",this_gene,".txt"))
# load gene expression and DEseq2 info
genestats <- fread(paste0(workingdir,"/gene_expression_and_DE_annotation_summary_wide.txt"))
genedetails <- fread(paste0(workingdir,"/gene_expression_and_DE_annotation_detail_long.txt"))
colnames(genedetails) <- c("gene","celltype","condition","meanCPM","gene_DESeq2_log2FC", "gene_DESeq2_padj", "gene_category")
if (nrow(genedetails) == 0) {
  genedetails <- data.frame(
    gene=character(),
    celltype=character(),
    condition=character(),
    meanCPM=numeric(),
    gene_DESeq2_log2FC=numeric(),
    gene_DESeq2_padj=numeric(),
    gene_category=character(),
    stringsAsFactors = FALSE
  )
}

# load genotype data
genotype <- fread(paste0(workingdir,"/genotype.txt"))
genotype.lean <- genotype %>% dplyr::select(c(ID,REF,ALT,CHROM,POS))
colnames(genotype.lean) <- c("snp","REF","ALT","snp_chr","snp_pos")

# load tested SNP ATACseq overlap info
atacoverlap <- fread(paste0(workingdir,"/snp_SNP_ATAC_overlap.txt")) %>%
  dplyr::select(-c("snp_chr","snp_start","snp_end")) %>%
  as.data.table()

atac_long <- melt(
  atacoverlap,
  id.vars = c("snp", "ATACpeak_chr", "ATACpeak_start", "ATACpeak_end","ATACpeak_name"),
  measure.vars = measure(value.name, celltype, condition,
    pattern = "^ATACpeak_(open|log2FoldChange|padj|Dynamic)_([^_]+)_([^_]+)$"
  )
)

colnames(atac_long) <- c("snp", "ATACpeak_chr", "ATACpeak_start", "ATACpeak_end", "ATACpeak_name", "celltype", "condition",
                         "ATACpeak_isopen", "ATACpeak_log2FC", "ATACpeak_padj")

if (nrow(atac_long) == 0) {
  atac_long <- data.frame(
    snp = character(),
    ATACpeak_chr = character(),
    ATACpeak_start = numeric(),
    ATACpeak_end = numeric(),
    ATACpeak_name = character(),
    celltype = character(),
    condition = character(),
    ATACpeak_isopen = character(),
    ATACpeak_log2FC = numeric(),
    ATACpeak_padj = numeric(),
    stringsAsFactors = FALSE
  )
}

## ---- 2. chromBPNet variant scores ----
scores <- fread(paste0(workingdir,"/chromBPNet_variant_scores.txt")) %>% as_tibble() %>%
  rename(snp = SNP) %>%
  dplyr::select(snp, celltype, condition,
                cbp_allele1 = allele1, cbp_allele2 = allele2,
                abs_logfc_x_active_allele_quantile.mean.pval,
                jsd_x_active_allele_quantile.mean.pval, 
                logfc.mean, jsd.mean, active_allele_quantile.mean) %>%    
  distinct(snp, celltype, condition, .keep_all = TRUE) %>%
  as.data.frame()

if (nrow(scores) == 0) {
  scores <- data.frame(
    snp = character(),
    celltype = character(),
    condition = character(),
    cbp_allele1 = character(),
    cbp_allele2 = character(),
    abs_logfc_x_active_allele_quantile.mean.pval = numeric(),
    jsd_x_active_allele_quantile.mean.pval = numeric(),
    logfc.mean = numeric(),
    jsd.mean = numeric(),
    active_allele_quantile.mean = numeric(),
    stringsAsFactors = FALSE
  )
}

## ---- 3. JOIN TABLES ----
df <- right_join(genotype.lean, modelstats, by="snp")

df1 <- left_join(df, scores, by=c("snp","celltype","condition"))

df2 <- left_join(df1, atac_long, by=c("snp", "celltype", "condition"))

df3 <- left_join(df2, genedetails, by=c("gene","celltype","condition")) %>%
  arrange(p)


## ---- 4. write output
fwrite(df3, file=paste0(workingdir,"/big_table.txt"), sep="\t", quote=F)



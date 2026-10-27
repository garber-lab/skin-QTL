#!/usr/bin/env Rscript
# fetch expression and cytokine DE info for a given gene

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
out_dir <- args[2]

# # for debug only
# this_gene <- "ERAP2"
# out_dir <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats/E/ERAP2"

# ------- annotate gene with expression and DE info and potential trans-acting role ------
path.Rscript.fetch_expressed_genes="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_expressed_genes.R"
source(path.Rscript.fetch_expressed_genes)
path.Rscript.fetch_DE_genes="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_DE_genes_log2FC1_padj0.05.R"
source(path.Rscript.fetch_DE_genes)
path.Rscript.fetch_induced_genes="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_induced_genes_log2FC1_padj0.05.R"
source(path.Rscript.fetch_induced_genes)
path.Rscript.fetch_mean_CPM_of_genes="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_mean_CPM_of_genes.R"
source(path.Rscript.fetch_mean_CPM_of_genes)
path.Rscript.fetch_DEseq_result_of_genes="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_DEseq_result_of_genes.R"
source(path.Rscript.fetch_DEseq_result_of_genes)
trans_acting_genes_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/literature/trans_acting_genes/compiled_trans_acting_candidate_genes_and_category.txt"
trans_acting_genes_table <- fread(trans_acting_genes_file) %>%
  set_colnames(c("gene","gene_category"))

gene_list <- this_gene
# --- high level summary ----
gene_exprs_annotation <- data.frame()
for (ct in c("FRB","KRT","MEL")) {
  for (cond in c("PBS", "IFNB", "IFNG", "TNF")) {
    
    expressed_genes <- unique(fetch_expressed_genes(ct, cond))
    
    if (cond == "PBS") {
      DE_genes <- character(0)
      induced_genes <- character(0)
    } else {
      DE_genes <- unique(fetch_DE_genes(ct, cond))
      induced_genes <- unique(fetch_induced_genes(ct, cond))
    }
    
    this_df <- data.frame(
      gene = gene_list,
      celltype = ct,
      condition = cond,
      is_expressed = gene_list %in% expressed_genes,
      is_DE = if (cond == "PBS") NA else gene_list %in% DE_genes,
      is_induced = if (cond == "PBS") NA else gene_list %in% induced_genes
    )
    
    gene_exprs_annotation <- dplyr::bind_rows(gene_exprs_annotation, this_df)
  }
}
gene_exprs_annotation <- left_join(gene_exprs_annotation, trans_acting_genes_table, by="gene")
fwrite(gene_exprs_annotation, file = paste0(out_dir,"/gene_expression_and_DE_annotation_summary_long.txt"), quote=F, sep="\t")

gene_exprs_annotation_wide <- gene_exprs_annotation %>%
  tidyr::pivot_wider(
    names_from = c(celltype, condition),
    values_from = c(is_expressed, is_DE, is_induced)
  )
fwrite(gene_exprs_annotation_wide, file = paste0(out_dir,"/gene_expression_and_DE_annotation_summary_wide.txt"), quote=F, sep="\t")

# --- detailed info ----
gene_detailed_annotation <- data.frame()
for (ct in c("FRB","KRT","MEL")) {
  for (cond in c("PBS", "IFNB", "IFNG", "TNF")) {
    
    # fetch mean CPM
    this_cpm <- fetch_mean_CPM_of_genes(ct, cond) %>%
      dplyr::filter(gene==this_gene) %>%
      dplyr::mutate(celltype = ct,
                    condition = cond) %>%
      dplyr::select(gene, celltype, condition, meanCPM)
    
    # fetch DESeq2 result
    if (cond=="PBS"){
      this_DESeq2_result <- data.frame(
        gene=this_gene,
        celltype=ct,
        condition=cond,
        log2FoldChange=NA,
        padj=NA)
    } else {
      this_DESeq2_result <- fetch_DEseq_result_of_genes(ct, cond) %>%
        dplyr::filter(gene==this_gene)
    }
    
    this_df <- left_join(this_cpm, this_DESeq2_result, by=c("gene","celltype","condition"))
    
    gene_detailed_annotation <- dplyr::bind_rows(gene_detailed_annotation, this_df)
    
    message(paste0(ct,"_",cond))
  }
}
gene_detailed_annotation <- left_join(gene_detailed_annotation, trans_acting_genes_table, by="gene")
fwrite(gene_detailed_annotation, file = paste0(out_dir,"/gene_expression_and_DE_annotation_detail_long.txt"), quote=F, sep="\t")

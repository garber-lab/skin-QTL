#!/usr/bin/env Rscript
# fetch all QTL modeling stats for a given gene across all celltypes, conditions, and QTLtype.

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
# out_dir <- 

# ------- compile all tested cis-QTL SNP table for all celltype x condition -----

this_gene_modeling_stats <- data.frame()

for (ct in c("FRB","KRT","MEL")) {
  
  # find gene chunk in this celltype
  dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
  if (ct %in% c("KRT", "FRB")) {
    this_gene_chunks <- fread(paste0(dir,"/data/gene_chunk_dict.txt")) %>%
      dplyr::filter(gene == this_gene) %>% pull(chunk)
  } else if (ct == "MEL" ) {
    this_gene_chunks <- fread(paste0(dir,"/data/gene_chunk_dict_old.txt")) %>%
      dplyr::filter(gene == this_gene) %>% pull(chunk)
  }
  this_gene_chunks <- sprintf("%03d", this_gene_chunks)
  
  # get the modeling stats for each condition and QTLtype
  for (cond in c("PBS", "IFNB", "IFNG", "TNF")) {
    for (QTLtype in c("eQTL", "reQTL")) {
      
      if (cond=="PBS" && QTLtype=="reQTL") {
        next
      } 
      message(ct,"_",cond,"_",QTLtype)
 
      # --- gather modeling stats across chunks ---
      this_modeling_stats <- data.frame()
      for (chunk in this_gene_chunks) {
        dir_model_result  <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct,"/results_QC/")
        fname <- paste0("modeling_stats_postQC_",ct,"_",cond,"_",QTLtype,"_",chunk,".txt")
        fpath <- file.path(dir_model_result, fname)
        if (!file.exists(fpath)) next
        chunk_tbl <- fread(fpath) %>%
          dplyr::filter(gene == this_gene) %>%
          dplyr::select(-qc_genotypes, -z) %>%
          dplyr::mutate(celltype = ct, condition = cond, QTLtype = QTLtype)
        this_modeling_stats <- rbind(this_modeling_stats, chunk_tbl)   # inside the loop
      }
      if (nrow(this_modeling_stats) == 0) next   # gene genuinely absent in this context
      
      # --- eigenMT cutoff ---
      eigenMT_summary <- fread(paste0(dir,"/eigenMT/results/",ct,"_",cond,"_",QTLtype,".eigenMT.txt"))
      sig_eigen <- eigenMT_summary[q_gene < 0.05]
      
      if (nrow(sig_eigen) == 0) {
        message("No eGenes for ", ct, " ", cond, " ", QTLtype)
        this_modeling_stats <- this_modeling_stats %>%
          dplyr::mutate(p_nominal_cutoff = NA_real_, is_egene = FALSE, is_sig_snp = FALSE)
      } else {
        ub <- min(c(eigenMT_summary[q_gene >= 0.05, p_gene_eigenMT], 1), na.rm = TRUE)
        p_thresh <- (max(sig_eigen$p_gene_eigenMT) + ub) / 2
        gene_annot <- eigenMT_summary[, .(gene,
                                          p_nominal_cutoff = p_thresh / Meff,
                                          is_egene = q_gene < 0.05)]
        this_modeling_stats <- this_modeling_stats %>%
          dplyr::left_join(gene_annot, by = "gene") %>%
          dplyr::mutate(is_egene   = dplyr::coalesce(is_egene, FALSE),   # gene absent from eigenMT -> FALSE
                        is_sig_snp = is_egene & (p < p_nominal_cutoff))
      }
      
      # --- single append point ---
      this_gene_modeling_stats <- rbind(this_gene_modeling_stats, this_modeling_stats)     
}}}

fwrite(this_gene_modeling_stats, 
       paste0(out_dir,"/compiled_modeling_stats_",this_gene,".txt"),
       quote = FALSE,
       col.names = TRUE,
       sep = "\t")
message(paste0("result written to ", out_dir,"/compiled_modeling_stats_",this_gene,".txt"))

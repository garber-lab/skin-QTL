#!/usr/bin/env Rscript
# BH correction WITHIN the a-priori gene family for a celltype x condition x QTLtype:
#   PBS                                       -> genes EXPRESSED in that celltype
#   cytokine (IFNG/IFNB/TNF), eQTL or reQTL   -> genes EXPRESSED AND INDUCED
# Restrict the aggregated eigenMT table to that family, run BH WITHIN it, rank by q,
# and write NEW files. The aggregated eigenMT table itself is left untouched.
#
# Usage:
#   Rscript step4_QTL_result_correction_BH.R <infile.eigenMT.txt> 
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("\nUsage:\n  Rscript step4_QTL_result_correction_BH.R <infile> \n\n")
  quit(save = "no", status = 1)
}
infile   <- args[[1]]

# toy example
# infile <- "FRB_IFNG_eQTL.eigenMT.txt"
basename <- unlist(strsplit(infile, ".eigenMT.txt", fixed = TRUE))[1]   # e.g. FRB_IFNG_eQTL
parts    <- unlist(strsplit(basename, "_", fixed = TRUE))
ct <- parts[1]; condition <- parts[2]; QTLtype <- parts[3]

infile_dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/", ct, "/eigenMT/results")
results_dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/results")
data_dir    <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/data"

read_list <- function(path) {
  if (!file.exists(path)) stop("missing gene list: ", path)
  g <- trimws(readLines(path)); unique(g[nzchar(g)])
}

## ---- a-priori BH family ----
if (condition == "PBS") {
  family_genes <- read_list(file.path(data_dir, paste0("gene_list_", ct, "_PBS_expressed.txt")))
  fam_desc <- "expressed"
} else {
  expr <- read_list(file.path(data_dir, paste0("gene_list_", ct, "_", condition, "_expressed.txt")))
  indu <- read_list(file.path(data_dir, paste0("gene_list_", ct, "_", condition, "_induced.txt")))
  family_genes <- intersect(expr, indu)
  fam_desc <- "expressed & induced"
}

## ---- read aggregated eigenMT (do NOT modify this file) ----
f <- fread(file.path(infile_dir, infile))
f[, p_gene_eigenMT := -expm1(Meff * log1p(-pmin))]   # stable Sidak; idempotent heal

## ---- restrict to family (tested = finite p_gene), then BH WITHIN it ----
f_fam <- f[gene %in% family_genes & is.finite(p_gene_eigenMT)]
f_fam[, q_gene := p.adjust(p_gene_eigenMT, method = "BH")]
setorder(f_fam, q_gene, pmin)
f_fdr05 <- f_fam[q_gene < 0.05]

message(sprintf("[%s %s %s] family (%s): %d genes | tested in family: %d | eGenes FDR<0.05: %d",
                ct, condition, QTLtype, fam_desc, length(family_genes), nrow(f_fam), nrow(f_fdr05)))


## ---- write NEW files; leave the eigenMT table intact ----
fwrite(f_fam,   file.path(results_dir, paste0(basename, "_gene_BH_table.txt")),    sep = "\t", quote = FALSE)
fwrite(f_fdr05, file.path(results_dir, paste0(basename, "_gene_fdr05_table.txt")), sep = "\t", quote = FALSE)
writeLines(f_fdr05$gene, file.path(results_dir, paste0(basename, "_gene_fdr05_genelist.txt")))

summary_file <- file.path(results_dir, "eGene_summary_by_context.txt")

cat(
  sprintf(
    "[%s %s %s] family (%s): %d genes | tested in family: %d | eGenes FDR<0.05: %d\n",
    ct, condition, QTLtype, fam_desc, length(family_genes), nrow(f_fam), nrow(f_fdr05)
  ),
  file = summary_file,
  append = TRUE
)
message("wrote ", basename, "_gene_BH_table.txt (+ _fdr05_table / _fdr05_genelist); eigenMT table left intact")

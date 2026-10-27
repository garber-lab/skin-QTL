#!/usr/bin/env Rscript
# =============================================================================
# step5_compile_QTL_output_of_union_eGenes.R
#
# Like step5_compile_QTL_output_of_eGenes.R, but instead of slicing each context
# to its OWN eGenes, it slices EVERY context to the GLOBAL UNION of eGenes -- the
# set of genes that are an eGene in ANY context (any celltype x any condition,
# for this QTL type). Each context's per-variant output is written for that whole
# union set (whatever subset was tested in that context).
#
# Why: the standard step5 dumps per-variant output only for a context's own
# FDR<0.05 eGenes. That blocks colocalization for any gene outside that context's
# eGene set. To run coloc ACROSS CELLTYPES and ACROSS CYTOKINES, every context
# needs the cis profile of any gene that is an eGene somewhere. The associations
# already exist in the post-QC chunks (the gene was tested); they were just never
# written out. This writes them.
#
# Inputs:
#   - eGene lists (all contexts): <ROOT>/results/<ct>_<cond>_<qtl>_gene_fdr05_genelist.txt
#   - per-chunk post-QC:          <ROOT>/<ct>/results_QC/modeling_stats_postQC_<ct>_<cond>_<qtl>_<chunk>.txt
# Output per context (ct, cond):
#   <ROOT>/results/QTL_output_of_all_tested_variants_of_union_of_eGenes_<ct>_<cond>_<qtl>.txt
#
# Usage:
#   Rscript step5_compile_QTL_output_of_union_eGenes.R [ROOT] [QTL]
#   defaults: ROOT=<cluster path>  QTL=eQTL
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models"
QTL  <- if (length(args) >= 2) args[2] else "eQTL"

results_dir <- file.path(ROOT, "results")
if (!dir.exists(results_dir)) stop("results dir not found: ", results_dir)
stamp <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))

## ---- discover ALL contexts for this QTL type, build the GLOBAL eGene union ----
gl <- list.files(results_dir, pattern = sprintf("_%s_gene_fdr05_genelist\\.txt$", QTL), full.names = FALSE)
if (length(gl) == 0) stop("No *_", QTL, "_gene_fdr05_genelist.txt in ", results_dir, " -- run step4 (BH) first.")
contexts <- sub(sprintf("_%s_gene_fdr05_genelist\\.txt$", QTL), "", gl)   # e.g. FRB_IFNB, MEL_TNF

read_egenes <- function(stem) {
  g <- unique(trimws(readLines(file.path(results_dir, sprintf("%s_%s_gene_fdr05_genelist.txt", stem, QTL)))))
  g[nzchar(g)]
}
global_union <- sort(unique(unlist(lapply(contexts, read_egenes))))
stamp("contexts (", QTL, "): ", length(contexts), " [", paste(contexts, collapse=", "), "]")
stamp("GLOBAL eGene union across all contexts: ", length(global_union), " genes")

## ---- for each context, slice its post-QC chunks to the GLOBAL union ----
for (stem in contexts) {
  parts <- strsplit(stem, "_", fixed = TRUE)[[1]]
  if (length(parts) != 2) { stamp("  [skip] cannot parse context '", stem, "' into ct_cond"); next }
  ct <- parts[1]; cond <- parts[2]

  qc_dir  <- file.path(ROOT, ct, "results_QC")
  pattern <- sprintf("^modeling_stats_postQC_%s_%s_%s_.*\\.txt$", ct, cond, QTL)
  qc_files <- list.files(qc_dir, pattern = pattern, full.names = TRUE)
  if (length(qc_files) == 0) { stamp("  [", ct, " ", cond, "] no post-QC chunks (", pattern, ") -- skip"); next }

  sliced <- vector("list", length(qc_files))
  for (i in seq_along(qc_files)) {
    dt <- fread(qc_files[i])
    if (nrow(dt) == 0 || !("gene" %in% names(dt))) next
    sliced[[i]] <- dt[gene %in% global_union]
  }
  out <- rbindlist(sliced, use.names = TRUE, fill = TRUE)
  if (nrow(out) > 0) {
    out[, `:=`(celltype = ct, condition = cond, QTLtype = QTL)]
    setcolorder(out, c("celltype","condition","QTLtype",
                       setdiff(names(out), c("celltype","condition","QTLtype"))))
  }
  outfile <- file.path(results_dir,
             sprintf("QTL_output_of_all_tested_variants_of_union_of_eGenes_%s_%s_%s.txt", ct, cond, QTL))
  fwrite(out, outfile, sep = "\t")
  stamp("  [", ct, " ", cond, "] union genes ", length(global_union),
        " -> genes tested here ", uniqueN(out$gene), ", variant-gene rows ", nrow(out),
        " -> ", basename(outfile))
}
stamp("Done.")

#!/usr/bin/env Rscript
# =============================================================================
# step5_compile_QTL_output_of_eGenes.R
#
# For each context (celltype x condition x QTLtype), gather ALL tested (QC-passed)
# variant-gene pairs that belong to a SIGNIFICANT eGene, by slicing the per-chunk
# post-QC modeling tables down to the significant-eGene set, and concatenating
# across chunks.
#
# Inputs per context (ct, cond, qtl):
#   - significant eGenes:  <ROOT>/results/<ct>_<cond>_<qtl>_gene_fdr05_genelist.txt   (from step4)
#   - per-chunk post-QC:   <ROOT>/<ct>/results_QC/modeling_stats_postQC_<ct>_<cond>_<qtl>_<chunk>.txt  (from step2)
# Output per context:
#   <ROOT>/results/QTL_output_of_all_tested_variants_of_eGenes_<ct>_<cond>_<qtl>.txt
#
# Contexts are auto-discovered from the *_gene_fdr05_genelist.txt files present in
# <ROOT>/results, so this works unchanged for any sibling pipeline.
#
# Usage:
#   Rscript step5_compile_QTL_output_of_eGenes.R [ROOT]
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB"

results_dir <- file.path(ROOT, "results")
if (!dir.exists(results_dir)) stop("results dir not found: ", results_dir)

## ---- discover contexts from the significant-eGene lists ----
genelist_files <- list.files(results_dir, pattern = "_gene_fdr05_genelist\\.txt$", full.names = FALSE)
if (length(genelist_files) == 0) {
  stop("No *_gene_fdr05_genelist.txt in ", results_dir, " -- run step4 (BH) first.")
}
contexts <- sub("_gene_fdr05_genelist\\.txt$", "", genelist_files)   # e.g. FRB_IFNG_eQTL
message(sprintf("[%s] Found %d context(s): %s",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S"), length(contexts),
                paste(contexts, collapse = ", ")))

for (ctx in contexts) {
  parts <- strsplit(ctx, "_", fixed = TRUE)[[1]]
  if (length(parts) != 3) {
    message("  [skip] cannot parse context '", ctx, "' into ct_cond_qtl"); next
  }
  ct <- parts[1]; cond <- parts[2]; qtl <- parts[3]

  ## sig eGenes for this context
  sig_genes <- unique(trimws(readLines(file.path(results_dir, paste0(ctx, "_gene_fdr05_genelist.txt")))))
  sig_genes <- sig_genes[nzchar(sig_genes)]
  if (length(sig_genes) == 0) {
    message(sprintf("  [%s] 0 significant eGenes -- nothing to compile", ctx)); next
  }

  ## per-chunk post-QC tables for this context
  qc_dir  <- file.path(ROOT, ct, "results_QC")
  pattern <- sprintf("^modeling_stats_postQC_%s_%s_%s_.*\\.txt$", ct, cond, qtl)
  qc_files <- list.files(qc_dir, pattern = pattern, full.names = TRUE)
  if (length(qc_files) == 0) {
    message(sprintf("  [%s] no post-QC chunk files in %s (pattern %s) -- skipping", ctx, qc_dir, pattern))
    next
  }

  ## read + slice to sig eGenes, chunk by chunk
  sliced <- vector("list", length(qc_files))
  for (i in seq_along(qc_files)) {
    dt <- fread(qc_files[i])
    if (nrow(dt) == 0 || !("gene" %in% names(dt))) next
    sliced[[i]] <- dt[gene %in% sig_genes]
  }
  out <- rbindlist(sliced, use.names = TRUE, fill = TRUE)

  ## stamp context (keeps the LONG-format grain: snp x celltype x condition x QTLtype)
  if (nrow(out) > 0) {
    out[, `:=`(celltype = ct, condition = cond, QTLtype = qtl)]
    setcolorder(out, c("celltype", "condition", "QTLtype",
                       setdiff(names(out), c("celltype", "condition", "QTLtype"))))
  }

  outfile <- file.path(results_dir,
                       paste0("QTL_output_of_all_tested_variants_of_eGenes_", ctx, ".txt"))
  fwrite(out, outfile, sep = "\t")
  message(sprintf("  [%s] %d sig eGenes -> %d variant-gene rows across %d chunks -> %s",
                  ctx, length(sig_genes), nrow(out), length(qc_files), basename(outfile)))
}

message(sprintf("[%s] Done.", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

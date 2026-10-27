#!/usr/bin/env Rscript
# =============================================================================
# step5_union_compile_QTL_output_of_union_eGenes.R
#
# Like step5_union_compile_QTL_output_of_eGenes.R, but slices each context's
# post-QC modeling tables to the UNION of eGenes across the two conditions being
# compared (default IFNG u IFNB, eQTL), rather than to that context's OWN eGenes.
#
# Why: the standard step5 dumps per-variant output only for a condition's own
# FDR<0.05 eGenes, so an IFNG-only eGene has no IFNB per-variant profile (and vice
# versa) -- which blocks colocalization for the IFNG-only / IFNB-only categories.
# Because the UNION run tests the SAME gene family in both conditions, the
# other-condition association already exists in the post-QC chunks; it was just
# never written out. This script writes it, giving every union eGene a full cis
# profile in BOTH conditions so coloc can run on all categories.
#
# For each celltype ct that has both <ct>_IFNG_<qtl> and <ct>_IFNB_<qtl>
# significant-eGene lists:
#     union_genes = eGenes(ct,IFNG,qtl)  U  eGenes(ct,IFNB,qtl)
# then for cond in {IFNG, IFNB}: slice per-chunk post-QC tables to union_genes.
#
# Inputs per (ct, cond, qtl):
#   - eGene lists:   <ROOT>/results_test_union_of_induced_genes/<ct>_<cond>_<qtl>_gene_fdr05_genelist.txt
#   - per-chunk QC:  <ROOT>/<ct>/results_QC/modeling_stats_postQC_<ct>_<cond>_<qtl>_<chunk>.txt
# Output per (ct, cond):
#   <ROOT>/results_test_union_of_induced_genes/
#          QTL_output_of_all_tested_variants_of_UNION_eGenes_<ct>_<cond>_<qtl>.txt
#
# Usage:
#   Rscript step5_union_compile_QTL_output_of_union_eGenes.R [ROOT] [COND_A] [COND_B] [QTL]
#   defaults: ROOT=<cluster path>  COND_A=IFNG  COND_B=IFNB  QTL=eQTL
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

args   <- commandArgs(trailingOnly = TRUE)
ROOT   <- if (length(args) >= 1) args[1] else "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB"
COND_A <- if (length(args) >= 2) args[2] else "IFNG"
COND_B <- if (length(args) >= 3) args[3] else "IFNB"
QTL    <- if (length(args) >= 4) args[4] else "eQTL"

results_dir <- file.path(ROOT, "results_test_union_of_induced_genes")
if (!dir.exists(results_dir)) stop("results dir not found: ", results_dir)
stamp <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))

read_egenes <- function(ct, cond) {
  f <- file.path(results_dir, sprintf("%s_%s_%s_gene_fdr05_genelist.txt", ct, cond, QTL))
  if (!file.exists(f)) return(character(0))
  g <- unique(trimws(readLines(f))); g[nzchar(g)]
}

## ---- discover celltypes that have BOTH condition eGene lists ----
gl <- list.files(results_dir, pattern = sprintf("_%s_gene_fdr05_genelist\\.txt$", QTL))
cts_A <- sub(sprintf("_%s_%s_gene_fdr05_genelist\\.txt$", COND_A, QTL), "",
             grep(sprintf("_%s_%s_", COND_A, QTL), gl, value = TRUE))
cts_B <- sub(sprintf("_%s_%s_gene_fdr05_genelist\\.txt$", COND_B, QTL), "",
             grep(sprintf("_%s_%s_", COND_B, QTL), gl, value = TRUE))
celltypes <- intersect(cts_A, cts_B)
if (length(celltypes) == 0)
  stop("No celltypes with both ", COND_A, " and ", COND_B, " ", QTL, " eGene lists in ", results_dir)
stamp("celltypes: ", paste(celltypes, collapse = ", "),
      "   (union over ", COND_A, " u ", COND_B, ", ", QTL, ")")

slice_context <- function(ct, cond, union_genes) {
  qc_dir  <- file.path(ROOT, ct, "results_QC")
  pattern <- sprintf("^modeling_stats_postQC_%s_%s_%s_.*\\.txt$", ct, cond, QTL)
  qc_files <- list.files(qc_dir, pattern = pattern, full.names = TRUE)
  if (length(qc_files) == 0) { stamp("  [", ct, " ", cond, "] no post-QC chunks (", pattern, ") -- skip"); return(invisible()) }
  sliced <- vector("list", length(qc_files))
  for (i in seq_along(qc_files)) {
    dt <- fread(qc_files[i])
    if (nrow(dt) == 0 || !("gene" %in% names(dt))) next
    sliced[[i]] <- dt[gene %in% union_genes]
  }
  out <- rbindlist(sliced, use.names = TRUE, fill = TRUE)
  if (nrow(out) > 0) {
    out[, `:=`(celltype = ct, condition = cond, QTLtype = QTL)]
    setcolorder(out, c("celltype","condition","QTLtype",
                       setdiff(names(out), c("celltype","condition","QTLtype"))))
  }
  outfile <- file.path(results_dir,
             sprintf("QTL_output_of_all_tested_variants_of_UNION_eGenes_%s_%s_%s.txt", ct, cond, QTL))
  fwrite(out, outfile, sep = "\t")
  stamp("  [", ct, " ", cond, "] union genes ", length(union_genes),
        " -> genes w/ rows ", uniqueN(out$gene), ", variant-gene rows ", nrow(out),
        " -> ", basename(outfile))
}

for (ct in celltypes) {
  union_genes <- union(read_egenes(ct, COND_A), read_egenes(ct, COND_B))
  stamp(ct, ": ", COND_A, " eGenes ", length(read_egenes(ct, COND_A)),
        " + ", COND_B, " eGenes ", length(read_egenes(ct, COND_B)),
        " -> UNION ", length(union_genes))
  slice_context(ct, COND_A, union_genes)
  slice_context(ct, COND_B, union_genes)
}
stamp("Done.")

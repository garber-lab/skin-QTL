#!/usr/bin/env Rscript
# =============================================================================
# step4_union_BH.R  --  UNION-FAMILY variant of step4_QTL_result_correction_BH.R
#
# EXPERIMENT ("test_union_of_induced_genes"): call IFNG and IFNB eGenes on a
# SHARED gene universe = genes expressed-&-induced in EITHER IFNG OR IFNB (the
# union of the two per-condition induced families), instead of gating each
# condition to its own induced set ("test_cytokine_induced_genes_only"). This lets
# a gene induced only by IFNG still be tested/called as an IFNB eGene (and vice
# versa) -- i.e. it reproduces the "off-diagonal" eGenes and makes IFNG-vs-IFNB a
# fully shared testing matrix.
#
# WHY only step4 changes: step1 models the union of IFNG/IFNB expressed-&-induced
# genes and fits ALL conditions per pair; step2 (QC) and step3 (eigenMT) run on
# every QC-passed gene with no induction gate. The per-condition induced family is
# applied ONLY here. So this reads the SAME concatenated eigenMT tables the gated
# run uses (rootdir/<ct>/eigenMT/results) and just widens the BH family.
#
# OUTPUTS -> rootdir/results_test_union_of_induced_genes/  (the gated results dir
# is left untouched).
#
# Usage:
#   Rscript step4_union_BH.R <infile.eigenMT.txt>      # e.g. FRB_IFNB_eQTL.eigenMT.txt
# =============================================================================
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("\nUsage:\n  Rscript step4_union_BH.R <infile.eigenMT.txt>\n\n")
  quit(save = "no", status = 1)
}
infile   <- args[[1]]
basename <- unlist(strsplit(infile, ".eigenMT.txt", fixed = TRUE))[1]   # e.g. FRB_IFNB_eQTL
parts    <- unlist(strsplit(basename, "_", fixed = TRUE))
ct <- parts[1]; condition <- parts[2]; QTLtype <- parts[3]

gw_root <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models"                          # reuse genome-wide gene lists
rootdir <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB" # this pipeline's outputs
RESULTS_SUBDIR <- "results_test_union_of_induced_genes"       # union output namespace
infile_dir  <- paste0(rootdir, "/", ct, "/eigenMT/results")   # SHARED, family-independent eigenMT tables
results_dir <- paste0(rootdir, "/", RESULTS_SUBDIR)
data_dir    <- paste0(gw_root, "/data")                       # expressed/induced gene lists reused from genome-wide run
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# conditions that share the union universe in this experiment
UNION_CONDS <- c("IFNG", "IFNB")

read_list <- function(path) {
  if (!file.exists(path)) stop("missing gene list: ", path)
  g <- trimws(readLines(path)); unique(g[nzchar(g)])
}
per_cond_family <- function(ct, cond) {
  # expressed-&-induced family for ONE cytokine condition (the gated definition)
  expr <- read_list(file.path(data_dir, paste0("gene_list_", ct, "_", cond, "_expressed.txt")))
  indu <- read_list(file.path(data_dir, paste0("gene_list_", ct, "_", cond, "_induced.txt")))
  intersect(expr, indu)
}

## ---- a-priori BH family ----
if (condition == "PBS") {
  family_genes <- read_list(file.path(data_dir, paste0("gene_list_", ct, "_PBS_expressed.txt")))
  fam_desc <- "expressed"
} else if (condition %in% UNION_CONDS) {
  # UNION: expressed-&-induced in EITHER IFNG or IFNB
  family_genes <- Reduce(union, lapply(UNION_CONDS, function(cc) per_cond_family(ct, cc)))
  fam_desc <- paste0("expressed & induced in ", paste(UNION_CONDS, collapse = " OR "), " (UNION)")
} else {
  # fallback (e.g. TNF): original per-condition gating
  family_genes <- per_cond_family(ct, condition)
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

## ---- write NEW files into the union namespace; leave the eigenMT table intact ----
fwrite(f_fam,   file.path(results_dir, paste0(basename, "_gene_BH_table.txt")),    sep = "\t", quote = FALSE)
fwrite(f_fdr05, file.path(results_dir, paste0(basename, "_gene_fdr05_table.txt")), sep = "\t", quote = FALSE)
writeLines(f_fdr05$gene, file.path(results_dir, paste0(basename, "_gene_fdr05_genelist.txt")))

summary_file <- file.path(results_dir, "eGene_summary_by_context.txt")
cat(
  sprintf("[%s %s %s] family (%s): %d genes | tested in family: %d | eGenes FDR<0.05: %d\n",
          ct, condition, QTLtype, fam_desc, length(family_genes), nrow(f_fam), nrow(f_fdr05)),
  file = summary_file, append = TRUE
)
message("wrote ", basename, " BH/fdr05 tables to ", RESULTS_SUBDIR, "; eigenMT table left intact")

#!/usr/bin/env Rscript
# =============================================================================
# step5_get_sig_QTLs_of_union_eGenes.R
#
# Emit the master list of variants that are SIGNIFICANT eQTLs OR reQTLs in ANY
# context (the union over all eGenes / all celltype x condition x QTLtype). This
# is the fine-mapping candidate set: for each eGene the significant cis variants
# are the ones ChromBPNet then prioritizes.
#
# Significance is back-calculated from the pipeline's eigenMT(Sidak)+BH scheme
# (mirrors chromBPNet_analysis/src/07_finemap_sig_eqtl.R):
#   per context : p_BH = max(p_gene_eigenMT among eGenes with q_gene < alpha)  # exact BH boundary
#   per eGene g : c_g  = 1 - (1 - p_BH)^(1/Meff_g)                             # invert Sidak
#   sig QTL     : nominal p < c_g   (within that context's eGenes)
# A variant is kept if it is significant in >=1 context; its contexts/QTL types
# are recorded. NO peak/accessibility filter (peak occupancy is a downstream
# annotation, not an inclusion gate).
#
# Inputs (per context <stem> = <ct>_<cond>_<qtl>, discovered from BH tables):
#   <ROOT>/results/<stem>_gene_BH_table.txt
#   <ROOT>/results/QTL_output_of_all_tested_variants_of_eGenes_<stem>.txt
#   <ROOT>/results/genotype_of_eGene_tested_variants.txt          (coords; optional)
# Output:
#   <ROOT>/results/list_of_variants_sig_eQTL_or_reQTL_of_union_eGenes.txt
#
# Usage: Rscript step5_get_sig_QTLs_of_union_eGenes.R [ROOT] [alpha]
#   defaults: ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models  alpha=0.05
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
setDTthreads(0)

args  <- commandArgs(trailingOnly = TRUE)
ROOT  <- if (length(args) >= 1) args[1] else "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models"
ALPHA <- if (length(args) >= 2) as.numeric(args[2]) else 0.05
results_dir <- file.path(ROOT, "results")
if (!dir.exists(results_dir)) stop("results dir not found: ", results_dir)
stamp <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
outfile <- file.path(results_dir, "list_of_variants_sig_eQTL_or_reQTL_of_union_eGenes.txt")

## ---- discover ALL contexts (eQTL + reQTL) from the BH tables ----------------
bh_files <- list.files(results_dir, pattern = "_gene_BH_table\\.txt$", full.names = FALSE)
if (length(bh_files) == 0) stop("No *_gene_BH_table.txt in ", results_dir, " -- run step4 (BH) first.")
stems <- sub("_gene_BH_table\\.txt$", "", bh_files)          # e.g. KRT_PBS_eQTL, FRB_IFNG_reQTL
stamp("contexts discovered: ", length(stems), " (", paste(stems, collapse=", "), ")")

## ---- per-context significant-QTL calling -----------------------------------
call_sig <- function(stem) {
  bh <- fread(file.path(results_dir, paste0(stem, "_gene_BH_table.txt")),
              select = c("gene", "Meff", "p_gene_eigenMT", "q_gene"))
  eg <- bh[q_gene < ALPHA & is.finite(p_gene_eigenMT)]
  if (nrow(eg) == 0) { stamp("  [", stem, "] 0 eGenes -- skip"); return(NULL) }
  p_BH <- max(eg$p_gene_eigenMT)                              # exact BH boundary
  eg[, c_g := -expm1(log1p(-p_BH) / Meff)]                    # 1-(1-p_BH)^(1/Meff), stable
  qf <- file.path(results_dir, paste0("QTL_output_of_all_tested_variants_of_eGenes_", stem, ".txt"))
  if (!file.exists(qf)) { stamp("  [", stem, "] MISSING ", basename(qf), " -- skip"); return(NULL) }
  q <- merge(fread(qf, select = c("snp", "gene", "p")), eg[, .(gene, c_g)], by = "gene")
  s <- unique(q[p < c_g, .(snp)])
  stamp("  [", stem, "] eGenes=", nrow(eg), "  p_BH=", signif(p_BH, 3),
        "  sig variants=", nrow(s))
  if (nrow(s) == 0) return(NULL)
  s[, `:=`(stem = stem, qtl = sub(".*_", "", stem), cond = tstrsplit(stem, "_")[[2]])]
  s[]
}
sig <- rbindlist(lapply(stems, call_sig), use.names = TRUE)
if (is.null(sig) || nrow(sig) == 0) stop("No significant QTLs found across contexts.")

## ---- collapse to one row per variant ---------------------------------------
per_var <- sig[, .(
  n_sig_contexts = uniqueN(stem),
  sig_contexts   = paste(sort(unique(stem)), collapse = ","),
  qtltypes       = paste(sort(unique(qtl)),  collapse = ","),
  scoreable_now  = any(cond %in% c("PBS", "IFNG"))   # ChromBPNet models exist for PBS+IFNG only
), by = snp]

## ---- attach coordinates (best-effort) --------------------------------------
gfile <- file.path(results_dir, "genotype_of_eGene_tested_variants.txt")
if (file.exists(gfile)) {
  geno <- unique(fread(gfile, select = c("CHROM", "POS", "ID", "REF", "ALT")), by = "ID")
  setnames(geno, c("chr", "pos", "snp", "ref", "alt"))
  per_var <- merge(per_var, geno, by = "snp", all.x = TRUE)
  per_var[, variant_id := ifelse(is.na(chr), snp, paste(snp, chr, pos, ref, alt, sep = "_"))]
  n_nocoord <- per_var[is.na(chr), .N]
  if (n_nocoord > 0) stamp("WARNING: ", n_nocoord, " variants lack coords (not in genotype file)")
} else {
  stamp("NOTE: ", basename(gfile), " not found -- writing rsID only, no coords")
  per_var[, `:=`(chr = NA_character_, pos = NA_integer_, ref = NA_character_,
                 alt = NA_character_, variant_id = snp)]
}

setcolorder(per_var, c("variant_id", "snp", "chr", "pos", "ref", "alt",
                       "qtltypes", "n_sig_contexts", "scoreable_now", "sig_contexts"))
setorder(per_var, chr, pos, snp)
fwrite(per_var, outfile, sep = "\t")

stamp("TOTAL sig eQTL/reQTL variants (union of eGenes): ", nrow(per_var),
      "  | scoreable now (PBS/IFNG models): ", per_var[scoreable_now == TRUE, .N])
stamp("wrote ", outfile)

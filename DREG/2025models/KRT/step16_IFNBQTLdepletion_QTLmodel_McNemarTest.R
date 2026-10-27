#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# McNemar's test for IFNB eQTL DEPLETION on matched donors.
#
# Design: the eQTL pipeline was re-run for two conditions (default IFNG vs IFNB)
# on the SAME ~54 donors present in BOTH conditions and the SAME union gene set.
# Each gene is therefore a PAIRED observation: eGene (BH q_gene < fdr) yes/no in
# condition A and in condition B. McNemar's test asks whether the two conditions
# have different marginal eGene rates, using only the DISCORDANT genes:
#
#                      B eGene        B not-eGene
#   A eGene              a                b        <- A-only (IFNG-only)
#   A not-eGene          c                d        <- B-only (IFNB-only)
#
# H0: b == c. Depletion of B (IFNB) => b >> c. Because the correction is matched
# (same donors, same genotypes => identical per-gene Meff across conditions; see
# eigenmt-meff-rank-fix), this comparison isolates the p-value/biology signal and
# removes donor-count as a power confound.
#
# Usage:
#   Rscript step16_IFNBQTLdepletion_QTLmodel_McNemarTest.R \
#           <ct> <condition_A> <condition_B> <QTLtype> [fdr] [results_dir]
# Defaults: FRB IFNG IFNB eQTL 0.05 <standard cluster results path>
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table) })

log <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")

## ---- args -----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
ct       <- if (length(args) >= 1) args[[1]] else "FRB"
condA    <- if (length(args) >= 2) args[[2]] else "IFNG"
condB    <- if (length(args) >= 3) args[[3]] else "IFNB"
QTLtype  <- if (length(args) >= 4) args[[4]] else "eQTL"
fdr      <- if (length(args) >= 5) as.numeric(args[[5]]) else 0.05
results_dir <- if (length(args) >= 6) args[[6]] else
  paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/", ct,
         "/IFNBeGene_depletion_issue/eigenMT/results")

fileA <- file.path(results_dir, paste0(ct, "_", condA, "_", QTLtype, ".eigenMT.txt"))
fileB <- file.path(results_dir, paste0(ct, "_", condB, "_", QTLtype, ".eigenMT.txt"))

log("McNemar depletion test: ", ct, " ", condA, " vs ", condB, " (", QTLtype,
    "), FDR<", fdr)
log("  A: ", fileA)
log("  B: ", fileB)

## ---- load per-gene BH-corrected tables ------------------------------------
# These are the post-BHcorrection.R aggregate tables (one row per gene, with q_gene).
read_gene_table <- function(f, cond) {
  if (!file.exists(f))
    stop("Missing eigenMT results file for ", cond, ": ", f,
         "\n  -> run step16_..._eigenMT_perGeneChunk.sh then _BHcorrection.sh first.")
  d <- fread(f)
  if (!"q_gene" %in% names(d))
    stop("No 'q_gene' column in ", f, " -- run BHcorrection.R (BH step) before McNemar.")
  # one row per gene (defensive: keep the strongest / smallest q per gene)
  d <- d[order(gene, q_gene)][, .SD[1L], by = gene]
  d[, .(gene, q_gene, p_gene_eigenMT)]
}

A <- read_gene_table(fileA, condA)
B <- read_gene_table(fileB, condB)
setnames(A, c("q_gene", "p_gene_eigenMT"), c("q_A", "p_A"))
setnames(B, c("q_gene", "p_gene_eigenMT"), c("q_B", "p_B"))

## ---- pair genes (full outer join on gene) ---------------------------------
m <- merge(A, B, by = "gene", all = TRUE)
m[, tested_A := !is.na(q_A)]
m[, tested_B := !is.na(q_B)]
# eGene call; %in% TRUE-style guard so NA (untested) never becomes an eGene.
m[, egene_A := !is.na(q_A) & q_A < fdr]
m[, egene_B := !is.na(q_B) & q_B < fdr]

## ---- universe breakdown ---------------------------------------------------
n_both  <- m[tested_A & tested_B, .N]
n_Aonly <- m[tested_A & !tested_B, .N]
n_Bonly <- m[!tested_A & tested_B, .N]
log("Gene universe: ", nrow(m), " total | tested in both = ", n_both,
    " | ", condA, "-only tested = ", n_Aonly, " | ", condB, "-only tested = ", n_Bonly)
if (n_Aonly + n_Bonly > 0)
  log("  (genes testable in only one condition cannot be paired; excluded from the",
      " PRIMARY test, included in the SENSITIVITY test as non-eGene.)")

## ---- helper: build 2x2 + run tests ----------------------------------------
run_mcnemar <- function(eA, eB, label) {
  a <- sum( eA &  eB); b <- sum( eA & !eB)
  c <- sum(!eA &  eB); d <- sum(!eA & !eB)
  M <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
              dimnames = list(A = c("eGene","not"), B = c("eGene","not")))
  cat("\n================ ", label, " ================\n", sep = "")
  cat("2x2 (rows =", condA, ", cols =", condB, "):\n"); print(M)
  cat(sprintf("  both eGene (a)      : %d\n", a))
  cat(sprintf("  %s-only eGene (b)   : %d\n", condA, b))
  cat(sprintf("  %s-only eGene (c)   : %d\n", condB, c))
  cat(sprintf("  neither (d)         : %d\n", d))
  cat(sprintf("  eGenes: %s = %d, %s = %d\n", condA, a + b, condB, a + c))

  if (b + c == 0) {
    cat("  No discordant genes (b+c=0) -> McNemar not applicable.\n")
    return(data.table(analysis = label, a, b, c, d,
                      n_eGene_A = a + b, n_eGene_B = a + c,
                      mcnemar_chisq = NA_real_, mcnemar_p = NA_real_,
                      exact_binom_p = NA_real_, direction = "no discordant pairs"))
  }
  mc_cc <- stats::mcnemar.test(M, correct = TRUE)                  # continuity-corrected chi-sq
  ex    <- stats::binom.test(b, b + c, p = 0.5)                    # exact (gold standard)
  dir   <- if (b > c) sprintf("%s DEPLETED (%s-only %d > %s-only %d)", condB, condA, b, condB, c)
           else if (c > b) sprintf("%s enriched (%s-only %d > %s-only %d)", condB, condB, c, condA, b)
           else "symmetric (b == c)"
  cat(sprintf("  discordant ratio b/c : %.2f\n", ifelse(c > 0, b / c, Inf)))
  cat(sprintf("  McNemar chi-sq (cc)  : %.3f, p = %.3g\n", mc_cc$statistic, mc_cc$p.value))
  cat(sprintf("  exact binomial p     : %.3g   [PRIMARY]\n", ex$p.value))
  cat(sprintf("  direction            : %s\n", dir))
  data.table(analysis = label, a, b, c, d, n_eGene_A = a + b, n_eGene_B = a + c,
             mcnemar_chisq = unname(mc_cc$statistic), mcnemar_p = mc_cc$p.value,
             exact_binom_p = ex$p.value, direction = dir)
}

## ---- PRIMARY: strict matched universe (genes tested in BOTH) ---------------
both <- m[tested_A & tested_B]
res_primary <- run_mcnemar(both$egene_A, both$egene_B,
                           "PRIMARY: genes tested in BOTH conditions")

## ---- SENSITIVITY: full union, one-condition-absent counted as non-eGene ----
# Upper bound on discordance: a gene called only in the condition where it was
# testable is treated as non-eGene in the other. Preempts the reviewer question
# "what about genes that dropped out of one condition?" (see ifnb-egene-deficit).
res_sens <- run_mcnemar(m$egene_A, m$egene_B,
                        "SENSITIVITY: full union (untested = non-eGene)")

## ---- write outputs --------------------------------------------------------
tag <- paste0(ct, "_", condA, "_vs_", condB, "_", QTLtype, "_fdr", fdr)
summary_out <- file.path(results_dir, paste0("mcnemar_summary_", tag, ".txt"))
pergene_out <- file.path(results_dir, paste0("mcnemar_pergene_", tag, ".txt"))

res <- rbindlist(list(res_primary, res_sens))
fwrite(res, summary_out, sep = "\t")

m[, category := fifelse(egene_A & egene_B, "both",
                 fifelse(egene_A & !egene_B, paste0(condA, "_only"),
                  fifelse(!egene_A & egene_B, paste0(condB, "_only"),
                   fifelse(tested_A & tested_B, "neither",
                    fifelse(tested_A, paste0(condA, "_only_tested"),
                            paste0(condB, "_only_tested"))))))]
fwrite(m[order(category, q_A, q_B)], pergene_out, sep = "\t")

log("Wrote summary : ", summary_out)
log("Wrote per-gene: ", pergene_out)
log("Done.")


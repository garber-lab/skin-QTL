#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# McNemar depletion test on COMMONLY-INDUCED genes (induced in BOTH conditions).
#
# Standalone companion to step16_..._McNemarTest.R (which runs the UNION of induced
# genes). "Commonly induced" = the INTERSECTION of the two per-condition induced
# lists, computed here internally (no pre-made file needed).
#
# Why a separate universe: the union is dominated by condition-specific induced
# genes (in FRB: 409 IFNB-only vs 139 IFNG-only induced), and a gene induced in
# only one condition has eQTL power mostly in that condition -- so the union eGene
# count is mechanically tilted toward the broader inducer. Restricting to commonly
# induced genes removes that induction-breadth confound: both conditions induced
# every gene here, so both had comparable expression/power.
#
# eGene = BH q_gene < fdr, RE-COMPUTED within the commonly-induced family. The
# eigenMT per-gene p (p_gene_eigenMT) is set-independent; only the BH family
# changes, so this is post-hoc on the existing *.eigenMT.txt tables (no rerun).
#
# Usage:
#   Rscript ...McNemarTest_commonlyInduced.R <ct> <condA> <condB> <QTLtype> [fdr] \
#           [results_dir] [inducedA_file] [inducedB_file]
# Defaults: FRB IFNG IFNB eQTL 0.05 + standard cluster paths.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table) })
log <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")

args    <- commandArgs(trailingOnly = TRUE)
ct      <- if (length(args) >= 1) args[[1]] else "FRB"
condA   <- if (length(args) >= 2) args[[2]] else "IFNG"
condB   <- if (length(args) >= 3) args[[3]] else "IFNB"
QTLtype <- if (length(args) >= 4) args[[4]] else "eQTL"
fdr     <- if (length(args) >= 5) as.numeric(args[[5]]) else 0.05
base    <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/", ct,
                  "/IFNBeGene_depletion_issue")
results_dir <- if (length(args) >= 6) args[[6]] else file.path(base, "eigenMT/results")
inducedA <- if (length(args) >= 7) args[[7]] else
  file.path(base, paste0("data/expressed_and_induced_genes_", ct, "_", condA, ".txt"))
inducedB <- if (length(args) >= 8) args[[8]] else
  file.path(base, paste0("data/expressed_and_induced_genes_", ct, "_", condB, ".txt"))

fileA <- file.path(results_dir, paste0(ct, "_", condA, "_", QTLtype, ".eigenMT.txt"))
fileB <- file.path(results_dir, paste0(ct, "_", condB, "_", QTLtype, ".eigenMT.txt"))
log("McNemar (commonly-induced): ", ct, " ", condA, " vs ", condB, " (", QTLtype, "), FDR<", fdr)

## ---- load per-gene eigenMT tables + induced lists -------------------------
read_gene_table <- function(f, cond) {
  if (!file.exists(f)) stop("Missing eigenMT results file for ", cond, ": ", f,
                            "\n  -> run eigenMT + BHcorrection first.")
  d <- fread(f)
  if (!"q_gene" %in% names(d)) stop("No 'q_gene' column in ", f, " -- run BHcorrection.R first.")
  d[order(gene, q_gene)][, .SD[1L], by = gene][, .(gene, q_gene, p_gene_eigenMT)]
}
A <- read_gene_table(fileA, condA); setnames(A, c("q_gene","p_gene_eigenMT"), c("q_A","p_A"))
B <- read_gene_table(fileB, condB); setnames(B, c("q_gene","p_gene_eigenMT"), c("q_B","p_B"))

rdlist <- function(f) { if (!file.exists(f)) stop("missing induced list: ", f)
                        unique(fread(f, header = FALSE)[[1]]) }
indA <- rdlist(inducedA); indB <- rdlist(inducedB)
commonly <- intersect(indA, indB)
log("induced: ", condA, "=", length(indA), ", ", condB, "=", length(indB),
    " -> commonly induced (intersection) = ", length(commonly), " genes")

## ---- restrict to commonly-induced + RE-RUN BH within that family ----------
m <- merge(A, B, by = "gene", all = TRUE)
m[, tested_A := !is.na(p_A)][, tested_B := !is.na(p_B)]
m <- m[gene %in% commonly]
log(nrow(m), " commonly-induced genes present in the eigenMT tables; re-running BH within them.")
m[tested_A == TRUE, q_A := p.adjust(p_A, method = "BH")]
m[tested_B == TRUE, q_B := p.adjust(p_B, method = "BH")]
m[, egene_A := !is.na(q_A) & q_A < fdr]
m[, egene_B := !is.na(q_B) & q_B < fdr]

n_both  <- m[tested_A & tested_B, .N]
log("tested in both = ", n_both, " | ", condA, "-only tested = ",
    m[tested_A & !tested_B, .N], " | ", condB, "-only tested = ", m[!tested_A & tested_B, .N])

## ---- 2x2 + tests ----------------------------------------------------------
run_mcnemar <- function(eA, eB, label) {
  a <- sum( eA &  eB); b <- sum( eA & !eB)
  c <- sum(!eA &  eB); d <- sum(!eA & !eB)
  M <- matrix(c(a, b, c, d), 2, 2, byrow = TRUE,
              dimnames = list(A = c("eGene","not"), B = c("eGene","not")))
  cat("\n================ ", label, " ================\n", sep = "")
  cat("2x2 (rows =", condA, ", cols =", condB, "):\n"); print(M)
  cat(sprintf("  both=%d  %s-only(b)=%d  %s-only(c)=%d  neither=%d\n", a, condA, b, condB, c, d))
  cat(sprintf("  eGenes: %s = %d, %s = %d\n", condA, a + b, condB, a + c))
  if (b + c == 0) { cat("  No discordant genes -> McNemar not applicable.\n")
    return(data.table(analysis = label, a, b, c, d, n_eGene_A = a+b, n_eGene_B = a+c,
                      mcnemar_p = NA_real_, exact_binom_p = NA_real_, direction = "no discordant pairs")) }
  mc <- stats::mcnemar.test(M, correct = TRUE)
  ex <- stats::binom.test(b, b + c, 0.5)
  dir <- if (b > c) sprintf("%s DEPLETED (%s-only %d > %s-only %d)", condB, condA, b, condB, c)
         else if (c > b) sprintf("%s enriched (%s-only %d > %s-only %d)", condB, condB, c, condA, b)
         else "symmetric"
  cat(sprintf("  McNemar chi-sq(cc) p = %.3g | exact binomial p = %.3g [PRIMARY]\n", mc$p.value, ex$p.value))
  cat(sprintf("  direction: %s\n", dir))
  data.table(analysis = label, a, b, c, d, n_eGene_A = a+b, n_eGene_B = a+c,
             mcnemar_p = mc$p.value, exact_binom_p = ex$p.value, direction = dir)
}
both <- m[tested_A & tested_B]
res  <- run_mcnemar(both$egene_A, both$egene_B, "COMMONLY-INDUCED (tested in both)")

## ---- write outputs --------------------------------------------------------
tag <- paste0(ct, "_", condA, "_vs_", condB, "_", QTLtype, "_fdr", fdr, "_commonlyInduced")
fwrite(res, file.path(results_dir, paste0("mcnemar_summary_", tag, ".txt")), sep = "\t")
m[, category := fifelse(egene_A & egene_B, "both",
                 fifelse(egene_A & !egene_B, paste0(condA, "_only"),
                  fifelse(!egene_A & egene_B, paste0(condB, "_only"), "neither")))]
fwrite(m[order(category, q_A, q_B)], file.path(results_dir, paste0("mcnemar_pergene_", tag, ".txt")), sep = "\t")
log("wrote mcnemar_summary_", tag, ".txt and mcnemar_pergene_", tag, ".txt")


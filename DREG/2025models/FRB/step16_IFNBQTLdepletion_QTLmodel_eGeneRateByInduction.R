#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# eGene DISCOVERY RATE by induction category.
#
# Normalizes eGene counts to the number of TESTED genes in each induction
# category, and tests whether the between-condition eGene-rate difference depends
# on induction category (i.e. is the apparent IFNB "enrichment" just an
# induction/power effect concentrated in IFNB-induced genes?).
#
# eGene = genome-wide BH q_gene < fdr -- ONE consistent definition (from
# BHcorrection.R) -- then STRATIFIED by category so rates are comparable across
# categories. (This differs from the induced-in-both McNemar, where the subset was
# the BH family; both are valid and answer different questions.)
#
# Categories (from the two induced-gene lists): condA-only | both | condB-only.
#
# Outputs: rate table (category x condition); per-category paired McNemar; and a
# condition x category interaction test (GLMM with a gene random effect, so the
# gene pairing is respected; falls back to GLM with a caveat if lme4 is absent).
#
# NOTE: rates control for gene-set COMPOSITION, not per-gene detection POWER. A
# residual power difference (e.g. higher IFNB expression/variance for a gene) can
# still inflate a rate. The permutation null remains the calibration arbiter.
#
# Usage:
#   Rscript ...eGeneRateByInduction.R <ct> <condA> <condB> <QTLtype> [fdr] \
#           [results_dir] [inducedA_file] [inducedB_file]
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

log("eGene rate by induction: ", ct, " ", condA, " vs ", condB, " (", QTLtype, "), FDR<", fdr)

## ---- load per-gene genome-wide BH tables + induced lists ------------------
rd <- function(f, cond) {
  if (!file.exists(f)) stop("missing ", cond, " eigenMT table: ", f)
  d <- fread(f)
  if (!"q_gene" %in% names(d)) stop("no q_gene column in ", f, " -- run BHcorrection.R first.")
  d[order(gene, q_gene)][, .SD[1L], by = gene][, .(gene, q_gene)]
}
A <- rd(file.path(results_dir, paste0(ct, "_", condA, "_", QTLtype, ".eigenMT.txt")), condA)
B <- rd(file.path(results_dir, paste0(ct, "_", condB, "_", QTLtype, ".eigenMT.txt")), condB)
setnames(A, "q_gene", "q_A"); setnames(B, "q_gene", "q_B")

rdlist <- function(f) { if (!file.exists(f)) stop("missing induced list: ", f)
                        unique(fread(f, header = FALSE)[[1]]) }
indA <- rdlist(inducedA); indB <- rdlist(inducedB)
log("induced genes: ", condA, "=", length(indA), ", ", condB, "=", length(indB),
    ", both=", length(intersect(indA, indB)))

## ---- assemble per-gene table + induction category -------------------------
m <- merge(A, B, by = "gene", all = TRUE)
m[, tested_A := !is.na(q_A)][, tested_B := !is.na(q_B)]
m[, eG_A := !is.na(q_A) & q_A < fdr][, eG_B := !is.na(q_B) & q_B < fdr]
m[, `:=`(inA = gene %in% indA, inB = gene %in% indB)]
m[, category := fifelse(inA & inB, "both",
                fifelse(inA & !inB, paste0(condA, "_only"),
                 fifelse(!inA & inB, paste0(condB, "_only"), "neither_induced")))]
cats <- c(paste0(condA, "_only"), "both", paste0(condB, "_only"))
m[, category := factor(category, levels = c(cats, "neither_induced"))]
if (m[category == "neither_induced", .N] > 0)
  log("note: ", m[category == "neither_induced", .N],
      " tested genes are in neither induced list (excluded from category analysis).")

## ---- RATE TABLE: category x condition -------------------------------------
rate1 <- function(cond) {
  te <- if (cond == condA) m$tested_A else m$tested_B
  eg <- if (cond == condA) m$eG_A else m$eG_B
  data.table(category = m$category, tested = te, eG = eg)[tested == TRUE][
    , .(condition = cond, n_tested = .N, n_eGene = sum(eG),
        eGene_rate = round(mean(eG), 4)), by = category]
}
rate_tab <- rbind(rate1(condA), rate1(condB))[order(category, condition)]
cat("\n===== eGene RATE by induction category (genome-wide FDR<", fdr, ") =====\n", sep = "")
print(rate_tab)

## ---- per-category paired McNemar (condA vs condB, denominator fixed) -------
mcn <- function(cat) {
  s <- m[category == cat & tested_A & tested_B]
  a <- sum(s$eG_A & s$eG_B); b <- sum(s$eG_A & !s$eG_B); c <- sum(!s$eG_A & s$eG_B)
  data.table(category = cat, n = nrow(s), both = a,
             A_only = b, B_only = c,
             exact_p = if (b + c > 0) binom.test(b, b + c, 0.5)$p.value else NA_real_)
}
cat("\n===== per-category paired McNemar (", condA, " vs ", condB,
    "; b=", condA, "-only eGene, c=", condB, "-only eGene) =====\n", sep = "")
print(rbindlist(lapply(cats, mcn)))

## ---- condition x category interaction -------------------------------------
long <- rbind(
  m[tested_A == TRUE, .(gene, category, condition = condA, eGene = as.integer(eG_A))],
  m[tested_B == TRUE, .(gene, category, condition = condB, eGene = as.integer(eG_B))])
long <- long[category %in% cats]
long[, category  := droplevels(factor(category, levels = cats))]
long[, condition := relevel(factor(condition), ref = condA)]

cat("\n===== condition x induction-category interaction =====\n")
cat("(does the ", condB, "-vs-", condA, " eGene-rate difference depend on induction category?)\n", sep = "")
if (requireNamespace("lme4", quietly = TRUE)) {
  suppressMessages({
    f0 <- lme4::glmer(eGene ~ condition + category + (1 | gene), long, family = binomial)
    f1 <- lme4::glmer(eGene ~ condition * category + (1 | gene), long, family = binomial)
  })
  an <- anova(f0, f1)
  cat(sprintf("  GLMM (gene random effect) interaction LRT: chisq=%.3f df=%d p=%.3g\n",
              an$Chisq[2], an$Df[2], an$`Pr(>Chisq)`[2]))
} else {
  g0 <- glm(eGene ~ condition + category, long, family = binomial)
  g1 <- glm(eGene ~ condition * category, long, family = binomial)
  an <- anova(g0, g1, test = "LRT")
  cat(sprintf("  GLM interaction LRT (pairing IGNORED, rough): p=%.3g  [install lme4 for the paired GLMM]\n",
              an$`Pr(>Chi)`[2]))
}

out <- file.path(results_dir, paste0("eGene_rate_by_induction_", ct, "_", condA,
                                     "_vs_", condB, "_", QTLtype, "_fdr", fdr, ".txt"))
fwrite(rate_tab, out, sep = "\t")
log("wrote rate table: ", out)


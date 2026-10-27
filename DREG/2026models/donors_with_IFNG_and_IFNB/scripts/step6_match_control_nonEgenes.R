#!/usr/bin/env Rscript
# =============================================================================
# step6_match_control_nonEgenes.R
#
# Build a matched NON-eGene control set for each context (celltype x condition x
# QTLtype), then compile the QTL output (post-QC tested variants) for that control
# set -- the control analogue of step5.
#
# Controls are non-eGenes matched to each eGene so that a non-significant control
# is a GENUINE null (comparable expression, induction, AND power) rather than an
# underpowered near-miss.
#
#   Distance (minimize)  : baseline expression = mean log2(CPM+1) over PBS samples;
#                          + induction = DESeq2 log2FoldChange (STIMULATED only).
#                          PBS -> 1-D (baseline only). VST is a switchable baseline.
#   Caliper (must match)  : lead-SNP min-genotype count (power proxy ~ MAF x n) and
#                          n_snps_used (cis-window testability). RELATIVE tolerance.
#   Pool                  : step4 BH family with q_gene > POOL_QMIN (a relaxed
#                          buffer above 0.05 so near-misses are excluded);
#                          for PBS the family is EXPRESSED, for cytokines INDUCED.
#
# Matching = greedy nearest-neighbor (hardest eGene first), 1:1 WITHOUT replacement,
# among candidates passing both calipers. If none pass -> relax caliper (flag
# caliper_relaxed); if the pool is exhausted -> reuse nearest (flag control_reused).
#
# A DIAGNOSTICS section at the end reports, per context, reuse/relax rates,
# candidate counts, and covariate balance (SMD), and prints concrete recommendations
# when a caliper looks too tight. Written to control_matching_diagnostics.txt.
#
# Outputs per context, in <ROOT>/results/ :
#   control_matching_<ct>_<cond>_<qtl>.txt
#   QTL_output_of_all_tested_variants_of_control_nonEgenes_<ct>_<cond>_<qtl>.txt
# Plus <ROOT>/results/control_matching_diagnostics.txt
#
# Usage:
#   Rscript step6_match_control_nonEgenes.R [ROOT] [BASELINE_SOURCE=CPM|VST]
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(magrittr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
ROOT            <- if (length(args) >= 1) args[1] else "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB"
BASELINE_SOURCE <- if (length(args) >= 2) toupper(args[2]) else "CPM"   # CPM (default) or VST

## ------------------ TUNABLE matching parameters ------------------
POOL_QMIN  <- 0.25   # controls must have q_gene > this (relaxed buffer above 0.05)
TOL_POWER  <- 0.20   # RELATIVE caliper on lead-SNP min-genotype count (+/-20%)
TOL_NSNPS  <- 0.50   # LOOSE relative caliper on n_snps_used (+/-50%)
# diagnostics warning thresholds
WARN_REUSE <- 0.05   # flag if >5% of eGenes reused a control (pool exhausted)
WARN_RELAX <- 0.20   # flag if >20% of eGenes needed the caliper relaxed
WARN_SMD   <- 0.25   # flag covariate imbalance if |SMD| exceeds this
WARN_ZEROC <- 0.20   # flag if >20% of eGenes had 0 in-caliper candidates

## ------------------ Paths ------------------
results_dir <- file.path(ROOT, "results")
CPM_FILE  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"
VST_FILE  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
META_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
FETCH_SCRIPT <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_DEseq_result_of_genes.R"

set.seed(42)
if (!dir.exists(results_dir)) stop("results dir not found: ", results_dir)
if (!file.exists(FETCH_SCRIPT)) stop("fetch_DEseq_result_of_genes.R not found: ", FETCH_SCRIPT)
source(FETCH_SCRIPT)

## ------------------ Baseline expression per celltype (PBS-only) ------------------
message(sprintf("[%s] Loading metadata + baseline (%s, PBS samples) ...",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S"), BASELINE_SOURCE))
meta <- fread(META_FILE)
stopifnot(all(c("sample","condition","celltype") %in% names(meta)))

baseline_for_ct <- function(ct) {
  pbs_samples <- unique(meta[celltype == ct & condition == "PBS", sample])
  if (length(pbs_samples) == 0) stop("No PBS samples for celltype ", ct)
  if (BASELINE_SOURCE == "VST") {
    m <- fread(VST_FILE); gcol <- "gene"
    samp_cols <- intersect(pbs_samples, setdiff(names(m), gcol))
    if (length(samp_cols) == 0) stop("No PBS VST columns for ", ct)
    bl <- rowMeans(as.matrix(m[, ..samp_cols]), na.rm = TRUE)
    out <- data.table(final_gene = m[[gcol]], baseline = bl)
  } else {
    m <- fread(CPM_FILE); gcol <- "final_gene"
    samp_cols <- intersect(pbs_samples, setdiff(names(m), c("gene","name","final_gene")))
    if (length(samp_cols) == 0) stop("No PBS CPM columns for ", ct)
    bl <- rowMeans(log2(as.matrix(m[, ..samp_cols]) + 1), na.rm = TRUE)
    out <- data.table(final_gene = m[[gcol]], baseline = bl)
  }
  out <- out[!is.na(final_gene) & nzchar(final_gene)]
  out[, .(baseline = mean(baseline, na.rm = TRUE)), by = final_gene]
}

## ------------------ Discover contexts ------------------
genelist_files <- list.files(results_dir, pattern = "_gene_fdr05_genelist\\.txt$", full.names = FALSE)
if (length(genelist_files) == 0) stop("No *_gene_fdr05_genelist.txt in ", results_dir, " -- run step4 first.")
contexts <- sub("_gene_fdr05_genelist\\.txt$", "", genelist_files)
message(sprintf("[%s] Found %d context(s). Params: POOL_QMIN=%.2f TOL_POWER=%.0f%% TOL_NSNPS=%.0f%%",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S"), length(contexts),
                POOL_QMIN, 100*TOL_POWER, 100*TOL_NSNPS))

baseline_cache <- list()
deg_cache      <- list()
match_all      <- list()   # accumulate per-context matches for diagnostics

# keep control within [target/(1+tol), target*(1+tol)]
within_rel <- function(cand, target, tol) cand >= target/(1+tol) & cand <= target*(1+tol)

for (ctx in contexts) {
  parts <- strsplit(ctx, "_", fixed = TRUE)[[1]]
  if (length(parts) != 3) { message("  [skip] cannot parse '", ctx, "'"); next }
  ct <- parts[1]; cond <- parts[2]; qtl <- parts[3]
  is_pbs <- (cond == "PBS")

  ## --- eGenes + non-eGene pool (from step4 BH table); q-buffer on the pool ---
  bh_file <- file.path(results_dir, paste0(ctx, "_gene_BH_table.txt"))
  if (!file.exists(bh_file)) { message("  [skip ", ctx, "] missing ", basename(bh_file)); next }
  bh <- fread(bh_file)
  if (!all(c("gene","q_gene","lead_snp","n_snps_used") %in% names(bh))) {
    message("  [skip ", ctx, "] BH table lacks gene/q_gene/lead_snp/n_snps_used"); next
  }
  eGenes   <- unique(bh[q_gene <  0.05,      gene])
  pool_all <- unique(bh[q_gene >  POOL_QMIN, gene])
  if (length(eGenes) == 0)   { message("  [", ctx, "] 0 eGenes -- skip"); next }
  if (length(pool_all) == 0) { message("  [", ctx, "] 0 non-eGene candidates at q>", POOL_QMIN, " -- skip"); next }

  ## --- read this context's post-QC once (for lead min-geno AND the compile) ---
  qc_dir <- file.path(ROOT, ct, "results_QC")
  qc_files <- list.files(qc_dir, pattern = sprintf("^modeling_stats_postQC_%s_%s_%s_.*\\.txt$", ct, cond, qtl),
                         full.names = TRUE)
  if (length(qc_files) == 0) { message("  [skip ", ctx, "] no post-QC files in ", qc_dir); next }
  qc_all <- rbindlist(lapply(qc_files, fread), use.names = TRUE, fill = TRUE)
  if (!all(c("gene","snp","n_geno0","n_geno1","n_geno2") %in% names(qc_all))) {
    message("  [skip ", ctx, "] post-QC lacks gene/snp/n_geno*"); next
  }
  qc_all[, min_geno := pmin(n_geno0, n_geno1, n_geno2, na.rm = TRUE)]

  ## --- covariates: baseline (+log2FC), n_snps_used, lead-SNP min-geno ---
  if (is.null(baseline_cache[[ct]])) baseline_cache[[ct]] <- baseline_for_ct(ct)
  base_dt <- baseline_cache[[ct]]

  # lead-SNP min-geno per gene: join BH lead_snp to qc_all
  lead_geno <- merge(bh[, .(gene, lead_snp, n_snps = n_snps_used)],
                     qc_all[, .(gene, snp, lead_min_geno = min_geno)],
                     by.x = c("gene","lead_snp"), by.y = c("gene","snp"), all.x = TRUE)

  cov <- merge(data.table(gene = union(eGenes, pool_all)), base_dt,
               by.x = "gene", by.y = "final_gene", all.x = TRUE)
  if (is_pbs) {
    cov[, log2FC := NA_real_]; nn_cols <- "baseline"
  } else {
    ck <- paste0(ct, "_", cond)
    if (is.null(deg_cache[[ck]])) {
      d <- as.data.table(fetch_DEseq_result_of_genes(ct, cond))
      d <- d[!is.na(log2FoldChange), .(gene, log2FC = log2FoldChange)]
      deg_cache[[ck]] <- d[, .(log2FC = mean(log2FC, na.rm = TRUE)), by = gene]
    }
    cov <- merge(cov, deg_cache[[ck]], by = "gene", all.x = TRUE)
    nn_cols <- c("baseline", "log2FC")
  }
  cov <- merge(cov, lead_geno[, .(gene, n_snps, lead_min_geno)], by = "gene", all.x = TRUE)
  cov[, is_eGene := gene %in% eGenes]

  # require all matching covariates present (NN + caliper)
  need <- c(nn_cols, "lead_min_geno", "n_snps")
  n_missing <- cov[!complete.cases(cov[, need, with = FALSE]), .N]
  cov <- cov[complete.cases(cov[, need, with = FALSE])]
  eG   <- cov[is_eGene == TRUE]
  pool <- cov[is_eGene == FALSE]
  if (nrow(eG) == 0 || nrow(pool) == 0) {
    message(sprintf("  [%s] after covariate join: %d eGenes, %d pool (%d dropped) -- skip",
                    ctx, nrow(eG), nrow(pool), n_missing)); next
  }

  ## --- standardize NN covariate(s); greedy caliper-constrained NN match ---
  allc <- rbind(eG, pool)
  mu  <- sapply(nn_cols, function(cc) mean(allc[[cc]]))
  sdv <- sapply(nn_cols, function(cc) { s <- sd(allc[[cc]]); if (!is.finite(s) || s == 0) 1 else s })
  zmat <- function(dt) sweep(sweep(as.matrix(dt[, nn_cols, with = FALSE]), 2, mu, "-"), 2, sdv, "/")
  Zpool <- zmat(pool); ZeG <- zmat(eG)
  pool_lg <- pool$lead_min_geno; pool_ns <- pool$n_snps
  pool_used <- rep(FALSE, nrow(pool))
  eG_order <- order(-sqrt(rowSums(ZeG^2)), eG$gene)   # hardest (most extreme) first

  res <- vector("list", nrow(eG))
  for (kk in seq_along(eG_order)) {
    k <- eG_order[kk]
    cal <- within_rel(pool_lg, eG$lead_min_geno[k], TOL_POWER) &
           within_rel(pool_ns, eG$n_snps[k],        TOL_NSNPS)
    n_cand <- sum(cal & !pool_used)                    # unused candidates within both calipers
    relaxed <- FALSE; reused <- FALSE
    avail <- which(cal & !pool_used)
    if (!length(avail)) { avail <- which(!pool_used); relaxed <- TRUE }          # relax caliper (still unused)
    if (!length(avail)) {                                                        # pool exhausted -> reuse
      avail <- which(cal); if (!length(avail)) avail <- seq_len(nrow(pool)); reused <- TRUE
    }
    d <- sqrt(rowSums(sweep(Zpool[avail, , drop = FALSE], 2, ZeG[k, ], "-")^2))
    j <- avail[which.min(d)]
    if (!reused) pool_used[j] <- TRUE
    res[[k]] <- data.table(
      celltype = ct, condition = cond, QTLtype = qtl,
      eGene = eG$gene[k], eGene_baseline = eG$baseline[k], eGene_log2FC = eG$log2FC[k],
      eGene_lead_min_geno = eG$lead_min_geno[k], eGene_n_snps = eG$n_snps[k],
      control_gene = pool$gene[j], control_baseline = pool$baseline[j], control_log2FC = pool$log2FC[j],
      control_lead_min_geno = pool$lead_min_geno[j], control_n_snps = pool$n_snps[j],
      dist_standardized = min(d), n_candidates = n_cand,
      caliper_relaxed = relaxed, control_reused = reused
    )
  }
  match_dt <- rbindlist(res)
  fwrite(match_dt, file.path(results_dir, paste0("control_matching_", ctx, ".txt")), sep = "\t")
  match_all[[ctx]] <- match_dt

  ## --- compile control QTL output (reuse the loaded post-QC) ---
  control_genes <- unique(match_dt$control_gene)
  ctrl_out <- qc_all[gene %in% control_genes]
  if (nrow(ctrl_out) > 0) {
    ctrl_out <- copy(ctrl_out); ctrl_out[, min_geno := NULL]
    ctrl_out[, `:=`(celltype = ct, condition = cond, QTLtype = qtl)]
    setcolorder(ctrl_out, c("celltype","condition","QTLtype",
                            setdiff(names(ctrl_out), c("celltype","condition","QTLtype"))))
  }
  fwrite(ctrl_out, file.path(results_dir,
         paste0("QTL_output_of_all_tested_variants_of_control_nonEgenes_", ctx, ".txt")), sep = "\t")

  message(sprintf("  [%s] eGenes=%d pool=%d | matched=%d reused=%d relaxed=%d | median_cand=%.0f | ctrl rows=%d (covariate-dropped=%d)",
                  ctx, length(eGenes), length(pool_all), nrow(match_dt),
                  sum(match_dt$control_reused), sum(match_dt$caliper_relaxed),
                  median(match_dt$n_candidates), nrow(ctrl_out), n_missing))
}

# =============================================================================
# DIAGNOSTICS -- is any caliper too tight? per-context metrics + recommendations
# =============================================================================
message("\n[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] === matching diagnostics ===")
if (length(match_all) == 0) { message("No contexts matched; nothing to diagnose."); quit(save = "no") }
M <- rbindlist(match_all)

# standardized mean difference between eGene and control covariate values
smd <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  s <- sqrt((var(a) + var(b)) / 2); if (!is.finite(s) || s == 0) return(NA_real_)
  (mean(a) - mean(b)) / s
}
diag <- M[, .(
  n_eGenes          = .N,
  reuse_rate        = mean(control_reused),
  relax_rate        = mean(caliper_relaxed & !control_reused),
  pct_zero_cand     = mean(n_candidates == 0),
  median_candidates = as.numeric(median(n_candidates)),
  smd_baseline      = smd(eGene_baseline,      control_baseline),
  smd_log2FC        = smd(eGene_log2FC,        control_log2FC),
  smd_lead_min_geno = smd(eGene_lead_min_geno, control_lead_min_geno),
  smd_n_snps        = smd(eGene_n_snps,        control_n_snps)
), by = .(celltype, condition, QTLtype)]

# per-context recommendation string
recommend <- function(r) {
  msg <- character(0)
  if (r$reuse_rate > WARN_REUSE)
    msg <- c(msg, sprintf("reuse=%.0f%% (pool exhausted): lower POOL_QMIN (e.g. %.2f) or widen calipers",
                          100*r$reuse_rate, max(0.05, POOL_QMIN/2)))
  if (r$relax_rate > WARN_RELAX)
    msg <- c(msg, sprintf("relax=%.0f%% (caliper too tight): raise TOL_POWER->%.0f%% and/or TOL_NSNPS->%.0f%%",
                          100*r$relax_rate, 100*min(1, TOL_POWER*1.5), 100*min(1.5, TOL_NSNPS*1.5)))
  if (r$pct_zero_cand > WARN_ZEROC)
    msg <- c(msg, sprintf("%.0f%% eGenes had 0 in-caliper candidates: widen calipers or enlarge pool", 100*r$pct_zero_cand))
  for (nm in c("smd_baseline","smd_log2FC","smd_lead_min_geno","smd_n_snps")) {
    v <- r[[nm]]
    if (is.finite(v) && abs(v) > WARN_SMD)
      msg <- c(msg, sprintf("imbalance %s SMD=%.2f (>%.2f): %s", sub("smd_","",nm), v, WARN_SMD,
                            if (grepl("min_geno|n_snps", nm)) "tighten that caliper" else "enlarge pool / lower POOL_QMIN"))
  }
  if (!length(msg)) "OK" else paste(msg, collapse = " | ")
}
diag[, recommendation := sapply(seq_len(.N), function(i) recommend(diag[i]))]

# print + write
print(diag[, .(celltype, condition, QTLtype, n_eGenes,
               reuse_rate = round(reuse_rate,3), relax_rate = round(relax_rate,3),
               median_candidates, pct_zero_cand = round(pct_zero_cand,3),
               smd_baseline = round(smd_baseline,2), smd_log2FC = round(smd_log2FC,2),
               smd_lead_min_geno = round(smd_lead_min_geno,2), smd_n_snps = round(smd_n_snps,2))])

diag_file <- file.path(results_dir, "control_matching_diagnostics.txt")
fwrite(diag, diag_file, sep = "\t")

# global recommendation
g_reuse <- mean(diag$reuse_rate); g_relax <- mean(diag$relax_rate)
message(sprintf("\nGLOBAL: mean reuse=%.1f%%  mean relax=%.1f%%  (params POOL_QMIN=%.2f TOL_POWER=%.0f%% TOL_NSNPS=%.0f%%)",
                100*g_reuse, 100*g_relax, POOL_QMIN, 100*TOL_POWER, 100*TOL_NSNPS))
if (g_reuse > WARN_REUSE || g_relax > WARN_RELAX) {
  message("RECOMMEND: calipers/pool look tight on average -> try TOL_POWER=", round(min(1, TOL_POWER*1.5),2),
          ", TOL_NSNPS=", round(min(1.5, TOL_NSNPS*1.5),2), ", POOL_QMIN=", max(0.05, POOL_QMIN/2), " and re-run.")
} else {
  message("RECOMMEND: calipers look adequate (low reuse/relax). Per-context notes in ", basename(diag_file), ".")
}
message("Wrote ", diag_file)
message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Done.")

#!/usr/bin/env Rscript
# =============================================================================
# opt_fit_qtl_subset.R  --  covariate-optimization fitter
#
# Same interaction-QTL model as step1_run_QTL_chunk.R
#   vst ~ geno + condition + geno:condition + <expr covariates> + <geno PCs> + (1|donor)
# per-condition eQTL via emtrends; reQTL (cond vs PBS) via contrasts,
# but parameterized so the covariate-optimization driver can swap in:
#   * ANY expression-covariate file (PEER / raw PC / condition-removed PC, any K)
#   * ANY number of genotype PCs (nGPC)
#   * a per-FOLD pair/genotype subset (instead of a genome-wide chunk)
#
# DIFFERENCE FROM step1 (the whole point):
#   step1 loaded a single 20-factor PEER file and took the first nPEER columns
#   AFTER dropping PEER1-4 (condition dummies) -- i.e. it subset a larger run.
#   Here the covariate file already contains exactly the K wanted factors/PCs,
#   so we use ALL of its non-`sample` columns. No slicing, no PEER1-4 skip.
#
# Output schema is IDENTICAL to step1 so steps 2-4 (QC / eigenMT / BH) run
# unchanged on the results.
#
# Usage:
#   Rscript opt_fit_qtl_subset.R <CELLTYPE> <COVARIATE_FILE> <nGPC> \
#           <PAIRS_FILE> <GENO_FILE> <OUTPUT_FILE> [n_workers]
# =============================================================================
Sys.setenv(TZ = "America/New_York")
suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(dplyr)
  library(lme4); library(lmerTest); library(emmeans); library(future.apply)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript opt_fit_qtl_subset.R <CELLTYPE> <COVARIATE_FILE> <nGPC> <PAIRS_FILE> <GENO_FILE> <OUTPUT_FILE> [n_workers]")
}
ct         <- args[1]
COV_FILE   <- args[2]                 # sample + covariate columns (any method/K)
nGPC       <- as.integer(args[3])
chunk_file <- args[4]                 # per-fold pairs file
GENO_FILE  <- args[5]                 # per-fold genotype file
output_file<- args[6]
n_workers  <- if (length(args) >= 7) as.integer(args[7]) else 5L

## ------------------ Fixed data paths (same as step1) ------------------ ##
VST_FILE      <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
META_FILE     <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
GENO_PCS_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/ancestry_PCA_rerun2026/output/genotype_PCs_for_modeling_imputed.txt"

message("====== opt QTL fit start ======")
message("CT=", ct, " | cov=", COV_FILE, " | nGPC=", nGPC)
message("pairs=", chunk_file, " | geno=", GENO_FILE, " | out=", output_file)

# multicore = fork: workers share memory (copy-on-write), so large globals like
# geno_dt are NOT serialized. Raise the globals ceiling so future's size *check*
# doesn't abort (the objects aren't actually copied under fork).
options(future.globals.maxSize = 16 * 1024^3)
plan(multicore, workers = n_workers); set.seed(1)
conds_all <- c("PBS","IFNG","IFNB","TNF")

## ---------- 1) Load pairs / VST / meta ---------- ##
pairs <- fread(chunk_file, header = TRUE)
colnames(pairs)[4]  <- "gene"
colnames(pairs)[11] <- "snp"
npairs <- nrow(pairs)
pairs$distance <- pairs$SNP_end - pairs$gene_end
message("Pairs in fold: ", npairs)

# Prefer the tiny per-sub-chunk VST subset written next to the pairs file
# (vst_chunk_<tag>.tsv) -- a few-gene fread instead of the 409MB whole VST, which
# is what keeps per-job memory tiny. Fall back to the full VST if absent.
vst_subset <- file.path(dirname(chunk_file),
                        sub("^pairs_chunk_", "vst_chunk_", basename(chunk_file)))
VST_TO_READ <- if (file.exists(vst_subset)) vst_subset else VST_FILE
message("VST source: ", if (file.exists(vst_subset)) "per-sub-chunk subset" else "full VST")
expr_wide <- fread(VST_TO_READ) %>%
  dplyr::filter(gene %in% pairs$gene) %>%
  dplyr::select(c("gene", contains(ct)))

meta <- fread(META_FILE) %>% dplyr::filter(celltype == ct)
stopifnot(all(c("sample","donor","condition") %in% names(meta)))
expr_samples <- intersect(colnames(expr_wide)[-1], meta$sample)
expr_wide <- as.data.table(expr_wide)[, c("gene", expr_samples), with = FALSE]
meta <- meta[sample %in% expr_samples]
meta$condition <- factor(as.character(meta$condition), levels = conds_all)

expr_long <- melt(as.data.table(expr_wide), id.vars = "gene",
                  variable.name = "sample", value.name = "vst")
expr_long <- merge(expr_long, meta, by = "sample", all.x = TRUE, all.y = FALSE)
rm(expr_wide); invisible(gc())

## ---------- 2) Load EXPRESSION COVARIATES (any method/K) ---------- ##
# The file already holds exactly the K columns we want. Use ALL of them.
covs <- fread(COV_FILE)
setnames(covs, 1, "sample")
cov_cols <- setdiff(colnames(covs), "sample")   # <- ALL of them; no slicing
message("Expression covariates (", length(cov_cols), "): ",
        if (length(cov_cols)) paste(cov_cols, collapse = ", ") else "(none, K=0)")
if (length(cov_cols) > 0) {
  expr_long <- merge(expr_long, covs, by = "sample", all.x = TRUE)
  expr_long <- expr_long[complete.cases(expr_long[, ..cov_cols]), ]
}

## ---------- 3) Load genotype PCs (first nGPC; PCs are ordered) ---------- ##
pc_use <- character(0)
if (nGPC > 0 && file.exists(GENO_PCS_FILE)) {
  gpcs <- fread(GENO_PCS_FILE)
  stopifnot("donor" %in% names(gpcs))
  pc_cols_all <- setdiff(colnames(gpcs), c("FID","donor","donor_num"))
  if (nGPC > length(pc_cols_all)) stop("Requested nGPC=", nGPC, " > available ", length(pc_cols_all))
  pc_use <- pc_cols_all[seq_len(nGPC)]
  expr_long <- merge(expr_long, gpcs[, c("donor", pc_use), with = FALSE], by = "donor", all.x = TRUE)
} else if (nGPC > 0) {
  warning("GENO_PCS_FILE not found; proceeding without genotype PCs.")
}

## ---------- 4) Genotype table + fast access ---------- ##
geno_dt <- fread(GENO_FILE)
fixed_cols <- c("CHROM","POS","ID","REF","ALT")
geno_sample_cols <- setdiff(colnames(geno_dt), fixed_cols)
if (length(geno_sample_cols) == 0) stop("No sample GT columns in GENOTYPE.")
geno_keep_samples <- intersect(geno_sample_cols, unique(expr_long$donor))
if (length(geno_keep_samples) == 0) stop("No overlapping donor IDs between GENOTYPE and expression.")
geno_dt <- geno_dt[, c(fixed_cols, geno_keep_samples), with = FALSE]
donors_with_geno_and_expr <- intersect(unique(expr_long$donor), geno_keep_samples)

pairs     <- as.data.frame(pairs)
expr_long <- as.data.frame(expr_long)
geno_dt   <- as.data.frame(geno_dt)
geno_donors <- setdiff(colnames(geno_dt), fixed_cols)

fetch_geno_for_snp <- function(snp_id, donor_vec) {
  row_idx <- match(snp_id, geno_dt$ID)
  if (is.na(row_idx)) return(rep(NA_real_, length(donor_vec)))
  g <- as.numeric(geno_dt[row_idx, geno_donors, drop = TRUE]); names(g) <- geno_donors
  g[donor_vec]
}

## ---------- 5) Model formula ---------- ##
cov_terms <- c(cov_cols, pc_use)
cov_part  <- if (length(cov_terms) > 0) paste("+", paste(cov_terms, collapse = " + ")) else ""
fml_str   <- paste0("vst ~ geno + condition + geno:condition ", cov_part, " + (1 | donor)")
message("Model formula: ", fml_str)

## ---------- 6) Per-pair fit (identical logic + output to step1) ---------- ##
fit_one <- function(snp, this_gene) {
  cols_needed <- c("sample","donor","condition","vst", cov_terms)
  df <- tryCatch({
    idx <- expr_long$gene == this_gene
    if (!any(idx)) return(NULL)
    expr_long[idx, cols_needed, drop = FALSE]
  }, error = function(e) NULL)
  if (is.null(df)) return(NULL)

  df$condition <- factor(as.character(df$condition))
  present_conds <- intersect(conds_all, unique(df$condition))
  df <- df %>% dplyr::filter(condition %in% present_conds)
  if ("PBS" %in% present_conds) df$condition <- stats::relevel(df$condition, ref = "PBS")
  else df$condition <- factor(df$condition, levels = present_conds)

  df$geno <- fetch_geno_for_snp(snp, df$donor)
  ug <- na.omit(unique(df$geno)); if (length(ug) < 2) return(NULL)
  qc_genotypes <- length(ug)
  tab_gc <- table(df$geno, df$condition)

  safe_get <- function(mat, rowname) {
    required <- c("PBS","IFNG","IFNB","TNF")
    if (!(rowname %in% rownames(mat))) return(setNames(rep(0, length(required)), required))
    row <- mat[rowname, ]; out <- numeric(length(required)); names(out) <- required
    existing <- intersect(required, names(row)); out[existing] <- row[existing]; out
  }
  qc_n_per_condition_by_genotype <- paste0(c(
    paste(safe_get(tab_gc,"0"), collapse="_"),
    paste(safe_get(tab_gc,"1"), collapse="_"),
    paste(safe_get(tab_gc,"2"), collapse="_")), collapse="_")

  fit <- tryCatch(
    lme4::lmer(as.formula(fml_str), data = df, REML = FALSE,
               control = lme4::lmerControl(optimizer = "bobyqa", calc.derivs = FALSE,
                 optCtrl = list(maxfun = 1e5), check.conv.singular = "ignore")),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  tr <- tryCatch(emmeans::emtrends(fit, ~ condition, var = "geno"), error = function(e) NULL)
  if (is.null(tr)) return(NULL)
  tr_summ <- as.data.frame(summary(tr, infer = c(TRUE, TRUE)))

  re_tbl <- NULL
  if ("PBS" %in% tr_summ$condition && length(unique(tr_summ$condition)) > 1) {
    re_con <- tryCatch(emmeans::contrast(tr, "trt.vs.ctrl", ref = "PBS"), error = function(e) NULL)
    if (!is.null(re_con)) { re_tbl <- as.data.frame(summary(re_con)); re_tbl$cond <- sub(" - PBS$", "", re_tbl$contrast) }
  }
  base_tbl <- as.data.frame(emmeans::emmeans(fit, specs = ~ condition | geno, at = list(geno = 0)))

  out <- data.frame(
    snp = snp, gene = this_gene, dist = pairs[which(pairs$snp==snp & pairs$gene==this_gene),"distance"],
    eQTL_baseline_PBS=NA_real_, eQTL_beta_PBS=NA_real_, eQTL_p_PBS=NA_real_, eQTL_se_PBS=NA_real_,
    eQTL_baseline_IFNG=NA_real_, eQTL_beta_IFNG=NA_real_, eQTL_p_IFNG=NA_real_, eQTL_se_IFNG=NA_real_,
    eQTL_baseline_IFNB=NA_real_, eQTL_beta_IFNB=NA_real_, eQTL_p_IFNB=NA_real_, eQTL_se_IFNB=NA_real_,
    eQTL_baseline_TNF=NA_real_, eQTL_beta_TNF=NA_real_, eQTL_p_TNF=NA_real_, eQTL_se_TNF=NA_real_,
    reQTL_dbeta_IFNG=NA_real_, reQTL_p_IFNG=NA_real_, reQTL_se_IFNG=NA_real_,
    reQTL_dbeta_IFNB=NA_real_, reQTL_p_IFNB=NA_real_, reQTL_se_IFNB=NA_real_,
    reQTL_dbeta_TNF=NA_real_, reQTL_p_TNF=NA_real_, reQTL_se_TNF=NA_real_,
    qc_genotypes=NA_real_, qc_n_per_condition_by_genotype=NA_real_
  )
  for (cond in intersect(conds_all, tr_summ$condition)) {
    rowc <- tr_summ[tr_summ$condition == cond, ]; row_base <- base_tbl[base_tbl$condition == cond, ]
    out[[paste0("eQTL_baseline_", cond)]] <- row_base$emmean
    out[[paste0("eQTL_beta_", cond)]] <- rowc$geno.trend
    out[[paste0("eQTL_p_", cond)]] <- rowc$p.value
    out[[paste0("eQTL_se_", cond)]] <- rowc$SE
  }
  if (!is.null(re_tbl)) for (cond in intersect(c("IFNG","IFNB","TNF"), re_tbl$cond)) {
    rowr <- re_tbl[re_tbl$cond == cond, ]
    out[[paste0("reQTL_dbeta_", cond)]] <- rowr$estimate
    out[[paste0("reQTL_p_", cond)]] <- rowr$p.value
    out[[paste0("reQTL_se_", cond)]] <- rowr$SE
  }
  out[["qc_genotypes"]] <- qc_genotypes
  out[["qc_n_per_condition_by_genotype"]] <- qc_n_per_condition_by_genotype
  out
}

## ---------- 7) Parallel over pairs ---------- ##
t_start <- Sys.time()
res_list <- future.apply::future_lapply(seq_len(npairs),
  function(i) fit_one(pairs$snp[i], pairs$gene[i]), future.seed = TRUE)
res <- data.table::rbindlist(res_list, fill = TRUE)
res$qc_n_donors_with_geno_and_expr <- length(donors_with_geno_and_expr)
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
data.table::fwrite(res, output_file, sep = "\t")
message("Wrote ", nrow(res), " rows to ", output_file)
message("====== opt QTL fit done in ",
        round(difftime(Sys.time(), t_start, units="mins"),2), " min ======")

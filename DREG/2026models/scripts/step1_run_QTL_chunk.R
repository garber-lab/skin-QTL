#!/usr/bin/env Rscript
# Mixed-effects eQTL / reQTL per SNP-gene pair (sample-level PEERs)
# Inputs (paths are CLI-args or defaults below):
#   1) chunk_file: TSV with columns: snp   gene
#   2) output_file: TSV to write results
# Environment files (edit paths or pass via args 3..6):
#   - VST.txt (wide): gene_id | <sample...>
#   - sample_metadata.txt: sample, donor, condition
#   - PEERs.txt: sample, PEER1..PEERk  (SAMPLE-LEVEL)
#   - genotype_pcs.txt: donor, PC1, PC2, ... (DONOR-LEVEL)
#   - GENOTYPE: TSV with columns: CHR, POS, START, END, ID, REF, ALT, <sample1>, <sample2>, ...
#               GTs in 0,1,2 format
#
# Model (one fit per pair):
#   cpm ~ geno * condition + PEERs + (1 | donor)
#   eQTL per condition via emtrends; reQTL (cond vs PBS) via contrasts
# REQUIREMENTS:
#   - No cross-condition donor matching (e.g. IFNB/IFNG) is applied here.
#     Each condition's model term is fit using whatever donors have data for
#     that condition; a donor missing IFNB, say, still contributes its PBS/
#     IFNG/TNF rows. The matched-donor restriction lives only in the separate
#     IFNB-depletion analysis, not in this general per-chunk QTL script.
#   - Donors without genotype data are dropped naturally: fetch_geno_for_snp()
#     returns NA for them, and lme4::lmer's default na.action removes any row
#     with an NA in a model term before fitting.
Sys.setenv(TZ = "America/New_York")
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(future.apply)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript run_QTL_chunk.R <CELLTYPE> <nExprPC> <nGPC> <chunk_id>")
}

ct <- args[1] # MEL, KRT, FRB
nExprPC <- as.integer(args[2]) # number of condition-removed expression PCs (K) to include
nGPC <- as.integer(args[3]) # number of genotype PCs to include
chunk_id  <- args[4]

chunk_file <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/data/chunks/pairs_chunk_",chunk_id,".tsv")
GENO_FILE <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/data/chunks/genotype_pairs_chunk_",chunk_id,".tsv") # columns: CHR POS ID REF ALT sample1 sample2 ...
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/",ct)
output_file <- paste0(dir,"/results/result_",chunk_id,".tsv")

## ------------------ Paths ------------------
VST_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
META_FILE     <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"   # columns: sample, donor, condition, etc
# Condition-removed expression PCs, one file per K, from generate_covariates.R.
# The K-specific file already holds exactly K columns (exprPC1..exprPCK) -- no
# condition dummies, no slicing needed downstream.
EXPR_PC_FILE  <- sprintf("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_opt/covariates/expression_pcs_condremoved/%s/exprPC_K%02d.tsv", ct, nExprPC)  # columns: sample, exprPC1, exprPC2, ...
GENO_PCS_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/ancestry_PCA_rerun2026/output/genotype_PCs_for_modeling_imputed.txt"      # columns: FID, donor, donor_num, PC1..PC10

message("====== QTL chunk start ======")
message("Chunk:  ", chunk_file)
message("Output: ", output_file)

## ---------- 0) Settings ----------
plan(multicore, workers = 5)
set.seed(1)

conds_all <- c("PBS","IFNG","IFNB","TNF")  # allowed condition labels

## ---------- 1) Load once ---------- ####
message("[1/7] Loading pairs from chunk ...")
pairs <- fread(chunk_file, header = TRUE)
colnames(pairs)[4] <- "gene"
colnames(pairs)[11] <- "snp"
npairs <- nrow(pairs)
pairs$distance = pairs$SNP_end - pairs$gene_end
message("Pairs in chunk: ", npairs)

message("[2/7] Loading VST ...")
# Prefer the tiny per-chunk VST subset (vst_chunk_<id>.tsv, written by
# step0.4_chunk_vst.R next to the pairs chunk) over the full 409MB VST -- keeps
# per-job memory low and skips the repeated big read. Falls back to full VST.
VST_CHUNK <- sub("pairs_chunk_", "vst_chunk_", chunk_file)
VST_TO_READ <- if (file.exists(VST_CHUNK)) VST_CHUNK else VST_FILE
message("      VST source: ", if (file.exists(VST_CHUNK)) "per-chunk subset" else "full VST")
expr_wide <- fread(VST_TO_READ) %>%
  dplyr::filter(gene %in% pairs$gene) %>%
  dplyr::select(c("gene", contains(ct)))

message("[3/7] Loading sample metadata ...")
meta <- fread(META_FILE) %>% dplyr::filter(celltype==ct)
stopifnot(all(c("sample","donor","condition") %in% names(meta)))

expr_samples <- intersect(colnames(expr_wide)[-1], meta$sample)
expr_wide <- as.data.table(expr_wide)[, c("gene", expr_samples), with = FALSE]
meta <- meta[sample %in% expr_samples]
meta$condition <- factor(as.character(meta$condition), levels = conds_all)

message("[4/7] Reshaping VST and merging metadata ...")
expr_long <- melt(as.data.table(expr_wide), id.vars = "gene",
                  variable.name = "sample", value.name = "vst")
expr_long <- merge(expr_long, meta, by = "sample", all.x = TRUE, all.y = FALSE)
rm(expr_wide); invisible(gc())

message("[5/7] Loading condition-removed expression PCs (K=", nExprPC, ") from ", EXPR_PC_FILE, " ...")
exprpcs <- fread(EXPR_PC_FILE)
colnames(exprpcs)[1] <- "sample"
stopifnot("sample" %in% names(exprpcs))
exprpc_cols <- setdiff(colnames(exprpcs), "sample")   # exactly K cols (exprPC1..exprPCK); use ALL
stopifnot(length(exprpc_cols) == nExprPC)

expr_long <- merge(expr_long, exprpcs, by = "sample", all.x = TRUE)
expr_long <- expr_long[complete.cases(expr_long[, ..exprpc_cols]), ]

message("[6/7] Loading genotype PCs ...")
if (file.exists(GENO_PCS_FILE)) {
  gpcs <- fread(GENO_PCS_FILE)
  stopifnot("donor" %in% names(gpcs))
  pc_cols <- setdiff(colnames(gpcs), "donor")
  expr_long <- merge(expr_long, gpcs, by = "donor", all.x = TRUE)
} else {
  pc_cols <- character(0)
  warning("GENO_PCS_FILE not found; proceeding without genotype PCs.")
}

message("[7/7] Loading GENOTYPE table ...")
geno_dt <- fread(GENO_FILE)
fixed_cols <- c("CHROM","POS","ID","REF","ALT")

geno_sample_cols <- setdiff(colnames(geno_dt), fixed_cols)
if (length(geno_sample_cols) == 0) stop("No sample GT columns found in GENOTYPE.")
# expr_long is not pre-filtered to any special donor set (no matched-donor
# restriction in this script) -- this intersect just trims GENOTYPE down to
# donors that actually have expression data in this chunk, for efficiency.
geno_keep_samples <- intersect(geno_sample_cols, unique(expr_long$donor))
if (length(geno_keep_samples) == 0) stop("No overlapping donor IDs between GENOTYPE and expression samples.")
geno_dt <- geno_dt[, c(fixed_cols, geno_keep_samples), with = FALSE]

# Donors with both expression and genotype data available in this chunk --
# reported for QC traceability, not used to filter anything further. Actual
# per-pair donor inclusion still depends on that donor having a non-NA
# genotype call for the specific SNP being tested (handled in fetch_geno_for_snp
# / lme4's na.action at fit time).
donors_with_geno_and_expr <- intersect(unique(expr_long$donor), geno_keep_samples)

pairs     <- as.data.frame(pairs)
expr_long <- as.data.frame(expr_long)
geno_dt   <- as.data.frame(geno_dt)

## ---------- 2) Prep fast genotype access ---------- ####
fixed_cols  <- c("CHROM","POS","ID","REF","ALT")
geno_donors <- setdiff(colnames(geno_dt), fixed_cols)  # donors in genotype table

fetch_geno_for_snp <- function(snp_id, donor_vec) {
  row_idx <- match(snp_id, geno_dt$ID)
  if (is.na(row_idx)) {
    return(rep(NA_real_, length(donor_vec)))
  }
  row <- geno_dt[row_idx, , drop = FALSE]
  g <- as.numeric(row[1, geno_donors, drop = TRUE])
  names(g) <- geno_donors
  g[donor_vec]
}

## ---------- 3) Build model formula ---------- ####
# expression-PC file already holds exactly K columns -> use ALL of them (no slicing,
# no condition-dummy skipping). Genotype PCs: take the first nGPC of PC1..PC10.
pc_cols <- setdiff(colnames(gpcs), c("FID","donor","donor_num"))
pc_use <- pc_cols[seq_len(nGPC)]
cov_terms <- c(exprpc_cols, pc_use)
cov_part <- if (length(cov_terms) > 0) paste("+", paste(cov_terms, collapse = " + ")) else ""
fml_str <- paste0("vst ~ geno + condition + geno:condition ", cov_part, " + (1 | donor)")
message("Model formula: ", fml_str)

## ---------- 4) Per-pair fit function ---------- ####
conds_all <- c("PBS","IFNG","IFNB","TNF")

fit_one <- function(snp, this_gene) {
  cols_needed <- c("sample","donor","condition","vst", exprpc_cols, pc_cols)
  df <- tryCatch({
    idx <- expr_long$gene == this_gene
    if (!any(idx)) return(NULL)
    expr_long[idx, cols_needed, drop = FALSE]
  }, error = function(e) NULL)

  df$condition <- factor(as.character(df$condition))
  present_conds <- intersect(conds_all, unique(df$condition))
  df <- df %>% dplyr::filter(condition %in% present_conds)

  if ("PBS" %in% present_conds) {
    df$condition <- stats::relevel(df$condition, ref = "PBS")
  } else {
    df$condition <- factor(df$condition, levels = present_conds)
  }

  df$geno <- fetch_geno_for_snp(snp, df$donor)

  ug <- na.omit(unique(df$geno))
  if (length(ug) < 2) return(NULL)
  qc_genotypes <- length(ug)
  tab_gc <- table(df$geno, df$condition)

  safe_get <- function(mat, rowname) {
    required <- c("PBS", "IFNG", "IFNB", "TNF")
    if (!(rowname %in% rownames(mat))) {
      return(setNames(rep(0, length(required)), required))
    }
    row <- mat[rowname, ]
    out <- numeric(length(required))
    names(out) <- required
    existing <- intersect(required, names(row))
    out[existing] <- row[existing]
    return(out)
  }
  qc_n_per_condition_genotype0 = paste(safe_get(tab_gc,"0"), collapse="_")
  qc_n_per_condition_genotype1 = paste(safe_get(tab_gc,"1"), collapse="_")
  qc_n_per_condition_genotype2 = paste(safe_get(tab_gc,"2"), collapse="_")
  qc_n_per_condition_by_genotype = paste0(c(qc_n_per_condition_genotype0,
                                            qc_n_per_condition_genotype1,
                                            qc_n_per_condition_genotype2), collapse="_")

  fml <- as.formula(fml_str)
  fit <- tryCatch(
    lme4::lmer(fml, data = df, REML = FALSE,
               control = lme4::lmerControl(
                 optimizer = "bobyqa",
                 calc.derivs = FALSE,
                 optCtrl = list(maxfun = 1e5),
                 check.conv.singular = "ignore"
               )),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  tr <- tryCatch(emmeans::emtrends(fit, ~ condition, var = "geno"),
                 error = function(e) NULL)
  if (is.null(tr)) return(NULL)
  tr_summ <- as.data.frame(summary(tr, infer = c(TRUE, TRUE)))

  re_tbl <- NULL
  if ("PBS" %in% tr_summ$condition && length(unique(tr_summ$condition)) > 1) {
    re_con <- tryCatch(emmeans::contrast(tr, "trt.vs.ctrl", ref = "PBS"),
                       error = function(e) NULL)
    if (!is.null(re_con)) {
      re_tbl <- as.data.frame(summary(re_con))
      re_tbl$cond <- sub(" - PBS$", "", re_tbl$contrast)
    }
  }
  base_tbl <- as.data.frame(
    emmeans::emmeans(
      fit,
      specs = ~ condition | geno,
      at    = list(geno = 0)
    )
  )

  out <- data.frame(
    snp = snp, gene = this_gene, dist = pairs[which(pairs$snp==snp & pairs$gene==this_gene),"distance"],
    eQTL_baseline_PBS = NA_real_, eQTL_beta_PBS = NA_real_, eQTL_p_PBS = NA_real_, eQTL_se_PBS = NA_real_,
    eQTL_baseline_IFNG = NA_real_, eQTL_beta_IFNG = NA_real_, eQTL_p_IFNG = NA_real_, eQTL_se_IFNG = NA_real_,
    eQTL_baseline_IFNB = NA_real_, eQTL_beta_IFNB = NA_real_, eQTL_p_IFNB = NA_real_, eQTL_se_IFNB = NA_real_,
    eQTL_baseline_TNF = NA_real_, eQTL_beta_TNF = NA_real_,  eQTL_p_TNF  = NA_real_, eQTL_se_TNF = NA_real_,
    reQTL_dbeta_IFNG = NA_real_, reQTL_p_IFNG = NA_real_, reQTL_se_IFNG = NA_real_,
    reQTL_dbeta_IFNB = NA_real_, reQTL_p_IFNB = NA_real_, reQTL_se_IFNB = NA_real_,
    reQTL_dbeta_TNF  = NA_real_, reQTL_p_TNF  = NA_real_, reQTL_se_TNF = NA_real_,
    qc_genotypes = NA_real_,
    qc_n_per_condition_by_genotype = NA_real_
  )

  for (cond in intersect(conds_all, tr_summ$condition)) {
    rowc <- tr_summ[tr_summ$condition == cond, ]
    row_base <- base_tbl[base_tbl$condition == cond, ]
    out[[paste0("eQTL_baseline_", cond)]] <- row_base$emmean
    out[[paste0("eQTL_beta_", cond)]] <- rowc$geno.trend
    out[[paste0("eQTL_p_",    cond)]] <- rowc$p.value
    out[[paste0("eQTL_se_",    cond)]] <- rowc$SE
  }

  if (!is.null(re_tbl)) {
    for (cond in intersect(c("IFNG","IFNB","TNF"), re_tbl$cond)) {
      rowr <- re_tbl[re_tbl$cond == cond, ]
      out[[paste0("reQTL_dbeta_", cond)]] <- rowr$estimate
      out[[paste0("reQTL_p_",     cond)]] <- rowr$p.value
      out[[paste0("reQTL_se_",     cond)]] <- rowr$SE
    }
  }

  out[["qc_genotypes"]] <- qc_genotypes
  out[["qc_n_per_condition_by_genotype"]] <- qc_n_per_condition_by_genotype

  out
}

## ---------- 5) Parallel over pairs ---------- ####
message("Fitting ", npairs, " pairs in parallel (", future::nbrOfWorkers(), " workers) ...")
t_start <- Sys.time()
res_list <- future.apply::future_lapply(
  seq_len(npairs),
  function(i) fit_one(pairs$snp[i], pairs$gene[i]),
  future.seed = TRUE
)

res <- data.table::rbindlist(res_list, fill = TRUE)

# Chunk-level count of donors with both expression and genotype data
# available (see donors_with_geno_and_expr above). This is a QC/traceability
# figure only -- it does not gate which donors enter any individual fit.
res$qc_n_donors_with_geno_and_expr <- length(donors_with_geno_and_expr)

data.table::fwrite(res, output_file, sep = "\t")
message("Wrote ", nrow(res), " rows to ", output_file)
t_end <- Sys.time()
message("====== QTL chunk done ======")
message("Time for fitting ", npairs, " pairs: ", round(difftime(t_end, t_start, units="mins"),2), "mins")

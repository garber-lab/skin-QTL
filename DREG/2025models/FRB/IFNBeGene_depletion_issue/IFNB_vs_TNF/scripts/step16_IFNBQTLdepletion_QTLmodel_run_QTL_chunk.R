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
#   - ONLY KEEP THE DONORS THAT HAVE DATA IN BOTH IFNB AND TNF (matched-donor design)
#   - ONLY KEEP DONORS THAT HAVE BEEN GENOTYPED
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
  stop("Usage: Rscript run_QTL_chunk.R <CELLTYPE> <nPEER> <nGPC> <chunk_id>")
}

ct <- args[1] # MEL, KRT, FRB
nPEER <- as.integer(args[2]) # number of PEER factors to include
nGPC <- as.integer(args[3]) # number of genotype PCs to include
chunk_id  <- args[4]

dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct,"/IFNBeGene_depletion_issue/IFNB_vs_TNF")
chunk_file <- paste0(dir,"/chunks/pairs_chunk_",chunk_id,".tsv")
GENO_FILE <- paste0(dir,"/chunks/genotype_pairs_chunk_",chunk_id,".tsv") # columns: CHR POS ID REF ALT sample1 sample2 ...
output_file <- paste0(dir,"/results/result_",chunk_id,".tsv")

# # toy example for debugging
# ct <- "FRB"
# nPEER=10
# nGPC=2
# chunk_id="096"
# dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct,"/IFNBeGene_depletion_issue/IFNB_vs_TNF")
# chunk_file <- paste0(dir,"/chunks/pairs_chunk_",chunk_id,".tsv")
# GENO_FILE <- paste0(dir,"/chunks/genotype_pairs_chunk_",chunk_id,".tsv") # columns: CHR POS ID REF ALT sample1 sample2 ...
# output_file <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/KRT/results/test_result_",chunk_id,".tsv")

## ------------------ Paths ------------------
# Expression + metadata
#CPM_FILE      <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"
VST_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
META_FILE     <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"   # columns: sample, donor, condition, etc
PEER_FILE     <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_factors/peer_factors_",ct,"_PBS-IFNG-IFNB-TNF.tsv")             # columns: sample, PEER1, PEER2, ...
GENO_PCS_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/genotype_PCs_for_modeling_07242025.txt"      # columns: donor, PC1, PC2, ...

message("====== QTL chunk start ======")
message("Chunk:  ", chunk_file)
message("Output: ", output_file)

## ---------- 0) Settings ----------
# Parallel workers (5 workers is pre-determined through benchmarking)
plan(multicore, workers = 5)
set.seed(1)

conds_all <- c("PBS","IFNG","IFNB","TNF")  # allowed condition labels

## ---------- 1) Load once ---------- ####
# Load pairs
message("[1/7] Loading pairs from chunk ...")
pairs <- fread(chunk_file, header = TRUE)
colnames(pairs)[4] <- "gene"
colnames(pairs)[11] <- "snp"
npairs <- nrow(pairs)
pairs$distance = pairs$SNP_end - pairs$gene_end
message("Pairs in chunk: ", npairs)

message("[2/7] Loading VST ...")
expr_wide <- fread(VST_FILE) %>%
  dplyr::filter(gene %in% pairs$gene) %>%
  dplyr::select(c("gene", contains(ct)))

message("[3/7] Loading sample metadata ...")
meta <- fread(META_FILE) %>% dplyr::filter(celltype==ct)
stopifnot(all(c("sample","donor","condition") %in% names(meta)))

## ---------- 3b) Restrict to qualifying donors ---------- ####
# Requirement: keep a donor ONLY IF (a) it has samples in BOTH IFNB and TNF
# for this celltype, AND (b) it has genotype data available.
message("[3b/7] Restricting to donors with IFNB+TNF data and genotype data ...")

donors_ifnb <- unique(meta$donor[meta$condition == "IFNB"])
donors_tnf <- unique(meta$donor[meta$condition == "TNF"])
donors_both_conditions <- intersect(donors_tnf, donors_ifnb)

# Peek at genotype chunk file header only (cheap) to get genotyped donor IDs
fixed_cols_geno <- c("CHROM","POS","ID","REF","ALT")
geno_header <- fread(GENO_FILE, nrows = 0)
donors_genotyped <- setdiff(colnames(geno_header), fixed_cols_geno)
if (length(donors_genotyped) == 0) stop("No donor genotype columns found in GENOTYPE chunk file header.")

donors_qualified <- intersect(donors_both_conditions, donors_genotyped)

message("  Donors with IFNB data:            ", length(donors_ifnb))
message("  Donors with TNF data:            ", length(donors_tnf))
message("  Donors with BOTH IFNB+TNF:       ", length(donors_both_conditions))
message("  Donors genotyped:                 ", length(donors_genotyped))
message("  Donors qualified (both+genotyped):", length(donors_qualified))
if (length(donors_qualified) == 0) {
  stop("No donors satisfy the both-condition + genotyped requirement for this chunk.")
}

# Restrict metadata to qualifying donors. NOTE: this keeps ALL conditions
# (PBS/IFNG/IFNB/TNF) for a qualifying donor, not just IFNB/TNF rows --
# PBS/IFNG data for these donors is still needed for the reQTL baseline
# and the IFNG specificity control.
meta <- meta[donor %in% donors_qualified]

# Keep VST columns that have metadata (for this celltype, qualifying donors only)
expr_samples <- intersect(colnames(expr_wide)[-1], meta$sample)
expr_wide <- as.data.table(expr_wide)[, c("gene", expr_samples), with = FALSE]
meta <- meta[sample %in% expr_samples]
meta$condition <- factor(as.character(meta$condition), levels = conds_all)

# melt VST to long (gene_id, sample, VST) and merge metadata
message("[4/7] Reshaping VST and merging metadata ...")
expr_long <- melt(as.data.table(expr_wide), id.vars = "gene",
                  variable.name = "sample", value.name = "vst")
expr_long <- merge(expr_long, meta, by = "sample", all.x = TRUE, all.y = FALSE)
rm(expr_wide); invisible(gc())

# PEER factors at sample level
message("[5/7] Loading PEER (sample-level) ...")
peers <- fread(PEER_FILE)
colnames(peers)[1] <- "sample"
stopifnot("sample" %in% names(peers))
peer_cols <- setdiff(colnames(peers), "sample")
stopifnot(length(peer_cols) > 0)

# Merge PEERs onto sample rows directly (sample-level PEER merge)
expr_long <- merge(expr_long, peers, by = "sample", all.x = TRUE)
# Drop rows without PEERs if any remain
expr_long <- expr_long[complete.cases(expr_long[, ..peer_cols]), ]

# Genotype PCs (per donor)
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

# Load genotype table
message("[7/7] Loading GENOTYPE table ...")
geno_dt <- fread(GENO_FILE)
fixed_cols <- c("CHROM","POS","ID","REF","ALT")

# Identify genotype sample columns
geno_sample_cols <- setdiff(colnames(geno_dt), fixed_cols)
if (length(geno_sample_cols) == 0) stop("No sample GT columns found in GENOTYPE.")

# Keep only genotype donors that appear in expression data
# (expr_long$donor is already restricted to qualifying donors from step 3b,
# so this is now a consistency check rather than the primary filter)
geno_keep_samples <- intersect(geno_sample_cols, unique(expr_long$donor))
if (length(geno_keep_samples) == 0) stop("No overlapping donor IDs between GENOTYPE and expression samples.")
geno_dt <- geno_dt[, c(fixed_cols, geno_keep_samples), with = FALSE]

pairs     <- as.data.frame(pairs)
expr_long <- as.data.frame(expr_long)
geno_dt   <- as.data.frame(geno_dt)

## ---------- 2) Prep fast genotype access ---------- ####
fixed_cols  <- c("CHROM","POS","ID","REF","ALT")
geno_donors <- setdiff(colnames(geno_dt), fixed_cols)  # donors in genotype table

fetch_geno_for_snp <- function(snp_id, donor_vec) {
  # find row index for this SNP
  row_idx <- match(snp_id, geno_dt$ID)
  if (is.na(row_idx)) {
    return(rep(NA_real_, length(donor_vec)))
  }
  row <- geno_dt[row_idx, , drop = FALSE]
  
  # numeric vector of genotypes for all donors in geno_dt
  g <- as.numeric(row[1, geno_donors, drop = TRUE])
  names(g) <- geno_donors
  
  # return genotypes in the order of donor_vec
  g[donor_vec]
}


## ---------- 3) Build model formula ---------- ####
peer_cols <- setdiff(colnames(peers), c("sample","PEER1","PEER2","PEER3","PEER4"))
pc_cols <- setdiff(colnames(gpcs), c("FID","donor","donor_num"))
# subset to number of factors to include
peer_use <- peer_cols[seq_len(nPEER)]
pc_use <- pc_cols[seq_len(nGPC)]
cov_terms <- c(peer_use, pc_use)
cov_part <- if (length(cov_terms) > 0) paste("+", paste(cov_terms, collapse = " + ")) else ""
# final formula string (random intercept by donor)
fml_str <- paste0("vst ~ geno + condition + geno:condition ", cov_part, " + (1 | donor)")
message("Model formula: ", fml_str)

## ---------- 4) Per-pair fit function ---------- ####
conds_all <- c("PBS","IFNG","IFNB","TNF")

fit_one <- function(snp, this_gene) {
  # Subset expression for this gene and keep columns we need
  cols_needed <- c("sample","donor","condition","vst", peer_cols, pc_cols)
  df <- tryCatch({
    idx <- expr_long$gene == this_gene
    if (!any(idx)) return(NULL)
    expr_long[idx, cols_needed, drop = FALSE]
  }, error = function(e) NULL)
  
  # Ensure condition factor (drop unused)
  df$condition <- factor(as.character(df$condition))
  present_conds <- intersect(conds_all, unique(df$condition))
  df <- df %>% dplyr::filter(condition %in% present_conds)
  
  # Relevel baseline to PBS when available
  if ("PBS" %in% present_conds) {
    df$condition <- stats::relevel(df$condition, ref = "PBS")
  } else {
    df$condition <- factor(df$condition, levels = present_conds) # no PBS — reQTL will become NA
  }
  
  # Attach donor genotypes for this SNP (numeric 0/1/2)
  df$geno <- fetch_geno_for_snp(snp, df$donor)
  
  # QC: Require at least 2 genotype groups and some minor allele support
  ug <- na.omit(unique(df$geno))
  if (length(ug) < 2) return(NULL)
  qc_genotypes <- length(ug)
  
  # QC: require at least 5 donors in each genotype group in each condition
  tab_gc <- table(df$geno, df$condition)
  
  safe_get <- function(mat, rowname) {
    required <- c("PBS", "IFNG", "IFNB", "TNF")
    if (!(rowname %in% rownames(mat))) {
      return(setNames(rep(0, length(required)), required))
    }
    # Extract row
    row <- mat[rowname, ]
    # Prepare output vector
    out <- numeric(length(required))
    names(out) <- required
    # Fill in existing values, leave missing as 0
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
  
  # Fit mixed model
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
  
  # eQTL: slope of geno within each present condition
  tr <- tryCatch(emmeans::emtrends(fit, ~ condition, var = "geno"),
                 error = function(e) NULL)
  if (is.null(tr)) return(NULL)
  tr_summ <- as.data.frame(summary(tr, infer = c(TRUE, TRUE)))  # columns: condition, geno.trend, SE, df, t.ratio, p.value
  
  # reQTL: delta slope vs PBS (only if PBS present)
  re_tbl <- NULL
  if ("PBS" %in% tr_summ$condition && length(unique(tr_summ$condition)) > 1) {
    re_con <- tryCatch(emmeans::contrast(tr, "trt.vs.ctrl", ref = "PBS"),
                       error = function(e) NULL)
    if (!is.null(re_con)) {
      re_tbl <- as.data.frame(summary(re_con))
      # contrasts like "IFNG - PBS" -> "IFNG"
      re_tbl$cond <- sub(" - PBS$", "", re_tbl$contrast)
    }
  }
  
  # baseline of each condition at genotype=0
  base_tbl <- as.data.frame(
    emmeans::emmeans(
      fit,
      specs = ~ condition | geno,
      at    = list(geno = 0)
    )
  )
  
  # Assemble one-row output
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
  
  # Fill eQTL per present condition
  for (cond in intersect(conds_all, tr_summ$condition)) {
    rowc <- tr_summ[tr_summ$condition == cond, ]
    row_base <- base_tbl[base_tbl$condition == cond, ]
    out[[paste0("eQTL_baseline_", cond)]] <- row_base$emmean
    out[[paste0("eQTL_beta_", cond)]] <- rowc$geno.trend
    out[[paste0("eQTL_p_",    cond)]] <- rowc$p.value
    out[[paste0("eQTL_se_",    cond)]] <- rowc$SE
  }
  
  # Fill reQTL deltas vs PBS
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
### ---------- 5) Parallel over pairs: benchmark timing with different n pairs ####
#TEST_N <- 10
#npairs <- TEST_N
#message("Fitting ", npairs, " pairs in parallel (", future::nbrOfWorkers(), " workers) ...")
#t_start <- Sys.time()
#res_list <- future.apply::future_lapply(
#  seq_len(npairs),
#  function(i) fit_one(pairs$snp[i], pairs$gene[i]),
#  future.seed = TRUE
#)
#t_end <- Sys.time()
#message("Time for fitting ", npairs, " pairs: ", round(difftime(t_end, t_start, units="secs"),2), " seconds")
#res <- data.table::rbindlist(res_list, fill = TRUE)

## ---------- 5) Parallel over pairs ---------- ####
message("Fitting ", npairs, " pairs in parallel (", future::nbrOfWorkers(), " workers) ...")
t_start <- Sys.time()
res_list <- future.apply::future_lapply(
  seq_len(npairs),
  function(i) fit_one(pairs$snp[i], pairs$gene[i]),
  future.seed = TRUE
)

res <- data.table::rbindlist(res_list, fill = TRUE)

# Record chunk-level donor qualification count for traceability. This is
# constant across all rows in the chunk since donor qualification happens
# once, upstream of per-pair modeling.
res$qc_n_donors_qualified <- length(donors_qualified)

data.table::fwrite(res, output_file, sep = "\t")
message("Wrote ", nrow(res), " rows to ", output_file)
t_end <- Sys.time()
message("====== QTL chunk done ======")
message("Time for fitting ", npairs, " pairs: ", round(difftime(t_end, t_start, units="mins"),2), "mins")

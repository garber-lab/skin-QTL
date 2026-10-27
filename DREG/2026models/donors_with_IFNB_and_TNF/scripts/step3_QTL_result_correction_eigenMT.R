#!/usr/bin/env Rscript
# script for controlling within-gene FWER of eQTL or reQTL modeling stats (pval) using eigenMT
# MATCHED-DONOR variant: the LD donor panel is restricted to the IFNB & TNF
# (+ genotype) subset (step0.0), matching the donor set the QTL was fit on.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(future.apply)
  library(magrittr)
  library(purrr)
})

#### Input ####
# ------------------------------#
# Main: parse CLI args ####
# ------------------------------#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("
Usage:
  Rscript step3_QTL_result_correction_eigenMT.R <ct> <condition> <QTLtype> <chunk_id>

Example:
  Rscript step3_QTL_result_correction_eigenMT.R FRB TNF eQTL 076
", "\n")
  quit(save = "no", status = 1)
}

ct             <- args[[1]]
condition      <- args[[2]]
QTLtype        <- args[[3]]
chunk_id       <- args[[4]] # 000 to 314
this_condition <- condition

rootdir    <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNB_and_TNF" # this pipeline's outputs + its own chunks
dir        <- paste0(rootdir, "/", ct)
chunks_dir <- paste0(rootdir, "/data/chunks")   # pairs/genotype chunks built by step0
pair_file       <- paste0(chunks_dir, "/pairs_chunk_", chunk_id, ".tsv")
geno_file       <- paste0(chunks_dir, "/genotype_pairs_chunk_", chunk_id, ".tsv")
modelstats_file <- paste0(dir, "/results_QC/modeling_stats_postQC_", ct, "_", condition, "_", QTLtype, "_", chunk_id, ".txt")
meta_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
out_prefix <- paste(ct, condition, QTLtype, chunk_id, sep = "_")
n_workers <- 2
base_seed <- 42

# matched-donor subset (IFNB & TNF expr + genotype), from step0.0
donor_subset_file <- paste0(rootdir, "/data/donor_subset_", ct, ".txt")
donor_subset <- unique(trimws(readLines(donor_subset_file)))
donor_subset <- donor_subset[nzchar(donor_subset)]

# --- eigenMT effective-test settings ---------------------------------------
# eigenMT's published default window is 200 SNPs, calibrated for cohorts with
# hundreds of samples. Our matched-donor list is small, so a single gene-wide
# SNP x SNP correlation matrix is rank-deficient (rank <= n_donors-1): Meff
# collapses to ~n_donors-1 for any gene with n_snps >> n_donors, regardless of
# true LD. We instead use the standard eigenMT sliding-window estimator with a
# window strictly below the donor count, so every within-window LD matrix is full
# rank and eigenvalues are genuine. window is capped at n_donors-1 inside
# run_one_gene_eigenMT.
eigenmt_window     <- 50L     # SNPs per non-overlapping window (< n_donors)
eigenmt_var_thresh <- 0.99    # eigenMT default: PCs explaining 99% of genotype variance

message("=== eigenMT (matched-donor: IFNB & TNF) for ", ct, " ", condition, " ", QTLtype, ", chunk: ", chunk_id, " ===")

#### ----Load data sets -------- ####
#### --- Load chunk files --------- ####
# load modeling stats for the QTL of interest of chunk
modelstats_subset <- fread(modelstats_file) %>%
  dplyr::mutate(key = paste0(gene,"_",snp))
table(modelstats_subset$gene)

# load SNP:gene pairs
pairs <- fread(pair_file, header=T)
table(pairs$gene_name)

pairs <- pairs %>%
  dplyr::mutate(key = paste0(gene_name,"_",SNP_ID)) %>%
  dplyr::filter(key %in% modelstats_subset$key)

# load genotype of all SNPs that passed QC in chunk
genotype_all <- fread(geno_file) %>%
  dplyr::filter(ID %in% unique(modelstats_subset$snp))

#### --- Load metadata to define LD donor set --- ####
# LD panel must be built from the SAME donor set used upstream in the QTL, i.e.
# the matched-donor subset intersected with the per-QTLtype donor requirement.
meta <- read_tsv(meta_file, show_col_types = FALSE) %>%
  dplyr::filter(celltype==ct) %>%
  dplyr::filter(donor %in% donor_subset) %>%          # matched-donor restriction
  dplyr::select(sample, donor, condition)

if (QTLtype == "reQTL") {
  donors_ld <- meta %>%
    dplyr::filter(condition %in% c("PBS", this_condition)) %>%
    dplyr::group_by(donor) %>%
    dplyr::filter(dplyr::n_distinct(condition) == 2) %>%   # paired: has both PBS and C
    dplyr::ungroup() %>%
    dplyr::pull(donor) %>% unique()
} else {  # eQTL
  donors_ld <- meta %>%
    dplyr::filter(condition == this_condition) %>%
    dplyr::pull(donor) %>% unique()
}
message("donors_ld (", QTLtype, " ", this_condition, "): ", length(donors_ld), " donors")

# ----------- eigenMT function for 1 gene ------ ####
# compute eigenMT stats for one gene. The two effective-test helpers are defined
# INSIDE this function on purpose: run_one_gene_eigenMT is dispatched to future
# workers by future_lapply, and nesting guarantees the helpers travel with it
# (no reliance on future/globals auto-exporting top-level helper functions).
run_one_gene_eigenMT <- function(this_gene,
                                 pairs, genotype_all, modelstats_subset, p_col="p",
                                 donors_ld,
                                 ct, condition, QTLtype,
                                 window = 50L, var_thresh = 0.99) {

  # --- effective-test helpers (nested so they ship to future workers) ---
  # eigenMT per-window count: number of leading eigenvalues (PCs) needed to explain
  # `var_thresh` of the genotype-correlation variance (= trace = m, #SNPs in window).
  find_num_eigs <- function(eigenvalues, variance, var_thresh) {
    eigenvalues <- sort(eigenvalues, decreasing = TRUE)
    running_sum <- 0; counter <- 0L
    target <- variance * var_thresh
    while (running_sum < target && counter < length(eigenvalues)) {
      counter <- counter + 1L
      running_sum <- running_sum + eigenvalues[counter]
    }
    counter
  }
  # Windowed eigenMT Meff (Davis et al. 2016; joed3/eigenMT): SNPs ordered by
  # position, split into NON-overlapping windows of `window` SNPs. window < n_donors
  # => each within-window LD matrix is full rank => genuine eigenvalues. Per-window
  # counts summed. Plain Pearson correlation (no Ledoit-Wolf shrinkage; windows are
  # full-rank by construction, and shrinkage was rejected for distorting Meff).
  meff_eigenMT_windowed <- function(G_snp_by_donor, pos, window = 50L, var_thresh = 0.99) {
    M <- nrow(G_snp_by_donor)
    if (M < 2L) return(as.numeric(M))
    ord <- order(pos)
    G   <- G_snp_by_donor[ord, , drop = FALSE]
    meff <- 0; start <- 1L
    while (start <= M) {
      stop <- min(start + window - 1L, M)
      Gw <- G[start:stop, , drop = FALSE]
      m  <- nrow(Gw)
      if (m == 1L) { meff <- meff + 1; start <- stop + 1L; next }
      R <- suppressWarnings(stats::cor(t(Gw), use = "pairwise.complete.obs"))
      R[!is.finite(R)] <- 0; diag(R) <- 1             # guard pairwise-complete NAs
      ev <- pmax(eigen(R, symmetric = TRUE, only.values = TRUE)$values, 0)
      meff <- meff + find_num_eigs(ev, m, var_thresh)
      start <- stop + 1L
    }
    meff
  }

  # SNP list for this gene from pairs
  snps <- pairs %>% dplyr::filter(gene_name == this_gene) %>% dplyr::pull(SNP_ID) %>% unique()

  # Pull SNP-level p-values for this gene
  ms_g <- modelstats_subset %>%
    dplyr::filter(gene == this_gene, snp %in% snps) %>%
    dplyr::select(snp, gene, dplyr::all_of(p_col))

  pvals <- suppressWarnings(as.numeric(ms_g[[p_col]]))
  names(pvals) <- ms_g$snp
  pvals <- pvals[is.finite(pvals) & !is.na(pvals)]

  # Build LD matrix from genotype_all for SNPs with p-values
  fixed_cols <- c("CHROM","POS","ID","REF","ALT")
  geno_sample_cols <- setdiff(colnames(genotype_all), fixed_cols)
  donors_use <- intersect(donors_ld, geno_sample_cols)

  snps_use <- intersect(names(pvals), genotype_all$ID)

  geno_sub <- genotype_all %>%
    dplyr::filter(ID %in% snps_use) %>%
    dplyr::select(ID, POS, dplyr::all_of(donors_use))

  # enforce SNP order (match snps_use)
  geno_sub <- geno_sub[match(snps_use, geno_sub$ID), , drop = FALSE]

  G_snp_by_donor <- as.matrix(dplyr::select(geno_sub, -ID, -POS))
  storage.mode(G_snp_by_donor) <- "numeric"

  # drop SNPs with zero variance across donors (causes NA correlations)
  sds <- apply(G_snp_by_donor, 1, stats::sd, na.rm = TRUE)
  keep <- which(is.finite(sds) & sds > 0)

  if (length(keep) == 0) {
    return(tibble(
      celltype = ct, condition = condition, QTLtype = QTLtype, gene = this_gene,
      n_snps_total = length(snps), n_snps_used = 0L, n_donors_ld = length(donors_use),
      pmin = NA_real_, lead_snp = NA_character_, Meff = NA_real_, p_gene_eigenMT = NA_real_,
      stage = "no_var_snps"
    ))
  }
  # Final SNP list after genotype+variance filtering
  snps_final <- geno_sub$ID[keep]
  pos_final  <- geno_sub$POS[keep]

  # Recompute pmin and lead_snp on FINAL SNP set
  pvals_final <- pvals[snps_final]
  pmin <- min(pvals_final, na.rm = TRUE)
  lead <- names(which.min(pvals_final))

  # Restrict genotype matrix to final SNP set
  G_snp_by_donor <- G_snp_by_donor[keep, , drop = FALSE]

  # If <2 SNPs remain, Meff ~ 1 and Sidak reduces to pmin
  if (nrow(G_snp_by_donor) < 2) {
    Meff <- 1
    p_gene <- -expm1(Meff * log1p(-pmin))   # stable Sidak: 1-(1-pmin)^Meff cancels to 0 for tiny pmin
    return(tibble(
      celltype = ct, condition = condition, QTLtype = QTLtype, gene = this_gene,
      n_snps_total = length(snps), n_snps_used = nrow(G_snp_by_donor), n_donors_ld = length(donors_use),
      pmin = pmin, lead_snp = lead, Meff = Meff, p_gene_eigenMT = p_gene,
      stage = "lt2_var_snps"
    ))
  }

  # Effective number of independent tests via windowed eigenMT.
  # Windows are capped below the donor count so each within-window LD matrix is
  # full rank -- this is the fix for the single-matrix rank-deficiency collapse.
  win  <- min(window, length(donors_use) - 1L)
  Meff <- meff_eigenMT_windowed(G_snp_by_donor, pos_final,
                                window = win, var_thresh = var_thresh)
  p_gene <- -expm1(Meff * log1p(-pmin))   # stable Sidak: 1-(1-pmin)^Meff cancels to 0 for tiny pmin

  tibble(
    celltype = ct,
    condition = condition,
    QTLtype = QTLtype,
    gene = this_gene,
    n_snps_total = length(snps),
    n_snps_used = nrow(G_snp_by_donor),
    n_donors_ld = length(donors_use),
    pmin = pmin,
    lead_snp = lead,
    Meff = Meff,
    p_gene_eigenMT = p_gene,
    stage = "ok"
  )
}
# ----------- MAIN ---------------------------- ####
# Main: Run eigenMT on all genes in chunk across 2 workers
# ------------------------------#
#### --- Process for every gene (chunk) --------- ####
start_time <- Sys.time()

# output folder and filename
out_dir <- paste0(dir, "/eigenMT/", condition,"/", QTLtype)
out_file <- file.path(out_dir, paste0(out_prefix, ".eigenMT.tsv"))

future::plan(future::multisession, workers = n_workers)

gene_lst <- unique(pairs$gene_name)

results_list <- future.apply::future_lapply(
  gene_lst,
  function(g) {
    run_one_gene_eigenMT(
      this_gene = g,
      pairs = pairs,
      genotype_all = genotype_all,
      modelstats_subset = modelstats_subset,
      p_col = "p",
      donors_ld = donors_ld,
      ct = ct,
      condition = condition,
      QTLtype = QTLtype,
      window = eigenmt_window,
      var_thresh = eigenmt_var_thresh
    )
  },
  future.seed = TRUE
)

summary_out <- dplyr::bind_rows(results_list)

end_time <- Sys.time()

cat("Start:", format(start_time), "\n")
cat("End:  ", format(end_time), "\n")
cat("Elapsed:", as.numeric(difftime(end_time, start_time, units = "secs")), "seconds\n")

readr::write_tsv(summary_out, out_file)
message("Wrote: ", out_file)

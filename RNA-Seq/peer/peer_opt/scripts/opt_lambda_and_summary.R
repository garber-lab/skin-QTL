#!/usr/bin/env Rscript
# =============================================================================
# opt_lambda_and_summary.R  --  score one (setting, fold, cell type)
#
# After the fitter + steps 2-4 have run for a given (setting, fold, cell type)
# under its own base-dir, this computes:
#   * genomic inflation lambda per condition x QTLtype (from step2 QC p-values)
#   * eGene count at FDR<0.05 per condition x QTLtype (from step4 fdr05 genelist)
# and APPENDS one row per (cell_type, fold, method, K, nGPC, condition, QTLtype)
# to a master summary TSV. Runtime and status are passed in by the driver.
#
# Discovery count alone is deliberately NOT used to pick a winner -- lambda is
# recorded alongside so a later step can balance discovery, calibration,
# stability (across folds), and parsimony.
#
# Base-dir layout (created by the driver; rootdir for steps 2-4):
#   <base-dir>/<CT>/results_QC/modeling_stats_postQC_<CT>_<cond>_<QTLtype>_<fold>.txt
#   <base-dir>/results/<CT>_<cond>_<QTLtype>_gene_fdr05_genelist.txt
#   <base-dir>/results/<CT>_<cond>_<QTLtype>_gene_BH_table.txt
#
# Usage:
#   Rscript opt_lambda_and_summary.R --base-dir DIR --cell-type CT --fold TAG \
#       --method M --K K --nGPC N --summary-file FILE \
#       [--runtime-sec S] [--status ok] [--n-genes-tested G] [--n-pairs P]
# =============================================================================
suppressPackageStartupMessages(library(data.table))

parse_args <- function(a) {
  o <- list(`base-dir`=NULL, `cell-type`=NULL, fold=NULL, method=NULL, K=NULL,
            nGPC=NULL, `summary-file`=NULL, `runtime-sec`=NA, status="ok",
            `n-genes-tested`=NA, `n-pairs`=NA, `seed-genes-file`=NULL)
  i <- 1; while (i <= length(a)) { k <- sub("^--","",a[[i]])
    if (!(k %in% names(o))) stop("Unknown option --", k); o[[k]] <- a[[i+1]]; i <- i+2 }; o
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
for (r in c("base-dir","cell-type","fold","method","K","nGPC","summary-file"))
  if (is.null(opt[[r]])) stop("--", r, " is required")
ct <- opt$`cell-type`; fold <- opt$fold

# genomic inflation factor from a vector of p-values (1 df)
lambda_gc <- function(p) {
  p <- suppressWarnings(as.numeric(p)); p <- p[is.finite(p) & p > 0 & p <= 1]
  if (length(p) < 10) return(NA_real_)
  stats::median(stats::qchisq(1 - p, df = 1)) / stats::qchisq(0.5, df = 1)
}

qc_dir  <- file.path(opt$`base-dir`, ct, "results_QC")
res_dir <- file.path(opt$`base-dir`, "results")

# seed eGenes for the recovery metric (sensitivity vs the genome-wide run)
seed_genes <- character(0)
if (!is.null(opt$`seed-genes-file`) && file.exists(opt$`seed-genes-file`)) {
  seed_genes <- unique(trimws(readLines(opt$`seed-genes-file`)))
  seed_genes <- seed_genes[nzchar(seed_genes)]
}

# discover which (cond, QTLtype) were run for this fold from the QC filenames
qc_files <- list.files(qc_dir,
  pattern = paste0("^modeling_stats_postQC_", ct, "_.*_", fold, "\\.txt$"), full.names = TRUE)
if (!length(qc_files)) {
  message("WARNING: no QC files for ", ct, " ", fold, " under ", qc_dir)
}

rows <- list()
for (f in qc_files) {
  # parse cond + QTLtype from: modeling_stats_postQC_<CT>_<cond>_<QTLtype>_<fold>.txt
  stem <- sub("\\.txt$", "", basename(f))
  stem <- sub(paste0("^modeling_stats_postQC_", ct, "_"), "", stem)
  stem <- sub(paste0("_", fold, "$"), "", stem)
  parts <- strsplit(stem, "_", fixed = TRUE)[[1]]
  cond <- parts[1]; qtltype <- parts[2]

  qc <- tryCatch(fread(f), error = function(e) NULL)
  lam <- if (!is.null(qc) && "p" %in% names(qc)) lambda_gc(qc$p) else NA_real_
  n_tested_qc <- if (!is.null(qc)) nrow(qc) else NA_integer_
  n_genes_qc  <- if (!is.null(qc) && "gene" %in% names(qc)) uniqueN(qc$gene) else NA_integer_

  # eGenes from step4 fdr05 genelist (preferred) or BH table (q_gene<0.05)
  gl <- file.path(res_dir, paste0(ct, "_", cond, "_", qtltype, "_gene_fdr05_genelist.txt"))
  bh <- file.path(res_dir, paste0(ct, "_", cond, "_", qtltype, "_gene_BH_table.txt"))
  n_egenes <- NA_integer_; n_genes_in_family <- NA_integer_; egene_set <- character(0)
  if (file.exists(gl)) {
    egene_set <- trimws(readLines(gl)); egene_set <- egene_set[nzchar(egene_set)]
    n_egenes <- length(egene_set)
  } else if (file.exists(bh)) {
    t <- tryCatch(fread(bh), error = function(e) NULL)
    if (!is.null(t) && "q_gene" %in% names(t)) {
      egene_set <- t$gene[t$q_gene < 0.05]; n_egenes <- length(egene_set); n_genes_in_family <- nrow(t)
    }
  }
  if (file.exists(bh) && is.na(n_genes_in_family)) {
    t <- tryCatch(fread(bh), error = function(e) NULL); if (!is.null(t)) n_genes_in_family <- nrow(t)
  }

  # recovery of the genome-wide seed eGenes (sensitivity): of seeded genes that
  # were actually tested here, how many are re-called at FDR<0.05.
  n_seed_tested <- NA_integer_; n_seed_recovered <- NA_integer_; recovery_rate <- NA_real_
  if (length(seed_genes) > 0 && !is.null(qc)) {
    seed_tested <- intersect(seed_genes, unique(qc$gene))
    n_seed_tested <- length(seed_tested)
    n_seed_recovered <- length(intersect(seed_tested, egene_set))
    recovery_rate <- if (n_seed_tested > 0) round(n_seed_recovered / n_seed_tested, 4) else NA_real_
  }

  rows[[length(rows)+1]] <- data.table(
    cell_type = ct, fold = fold, method = opt$method,
    K = as.integer(opt$K), nGPC = as.integer(opt$nGPC),
    condition = cond, QTLtype = qtltype,
    n_genes_tested   = as.integer(opt$`n-genes-tested`),
    n_pairs_tested   = as.integer(opt$`n-pairs`),
    n_pairs_postQC   = n_tested_qc,
    n_genes_postQC   = n_genes_qc,
    n_genes_in_family= n_genes_in_family,
    n_eGenes_fdr05   = n_egenes,
    n_seed_tested    = n_seed_tested,
    n_seed_recovered = n_seed_recovered,
    recovery_rate    = recovery_rate,
    lambda_gc        = round(lam, 4),
    runtime_sec      = suppressWarnings(as.numeric(opt$`runtime-sec`)),
    status           = opt$status,
    timestamp        = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}

if (!length(rows)) {
  rows[[1]] <- data.table(
    cell_type=ct, fold=fold, method=opt$method, K=as.integer(opt$K),
    nGPC=as.integer(opt$nGPC), condition=NA, QTLtype=NA,
    n_genes_tested=as.integer(opt$`n-genes-tested`), n_pairs_tested=as.integer(opt$`n-pairs`),
    n_pairs_postQC=NA, n_genes_postQC=NA, n_genes_in_family=NA, n_eGenes_fdr05=NA,
    n_seed_tested=NA, n_seed_recovered=NA, recovery_rate=NA,
    lambda_gc=NA, runtime_sec=suppressWarnings(as.numeric(opt$`runtime-sec`)),
    status=paste0(opt$status, ";no_qc_files"), timestamp=format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}
out <- rbindlist(rows)

sf <- opt$`summary-file`
dir.create(dirname(sf), showWarnings = FALSE, recursive = TRUE)
append <- file.exists(sf)
fwrite(out, sf, sep = "\t", append = append, col.names = !append)
message("Appended ", nrow(out), " row(s) to ", sf)
print(out[, .(condition, QTLtype, n_eGenes_fdr05, n_seed_recovered, n_seed_tested, recovery_rate, lambda_gc)])

#!/usr/bin/env Rscript
# =============================================================================
# generate_covariates.R
#
# Generate expression-covariate matrices (PEER factors and/or expression PCs)
# for eQTL / interaction-QTL covariate optimization.
#
# KEY DESIGN POINTS (see README_covariates.md for the full rationale):
#   * PEER factors are fit SEPARATELY for each requested K. If you want 5
#     factors you fit a PEER model with Nk = 5; you never subset a larger run.
#     PEER factors are not ordered/orthogonal like PCs, so subsetting is wrong.
#   * Expression PCs ARE ordered/orthogonal, so they are computed once up to
#     max(K) and subset for each smaller K.
#   * Everything is done per cell type, pooling all conditions within a cell
#     type (matching the QTL model, which is fit across conditions with
#     condition + genotype:condition terms).
#   * Condition is passed to PEER as a fixed-effect covariate (model.matrix
#     ~0+condition) so PEER does not absorb the stimulation axis you model
#     explicitly. PEER_getX() then returns [condition covariates | K factors];
#     we STRIP the leading condition columns and write ONLY the K learned
#     factors as PEER1..PEERK. (This is why the old file needed PEER1-4 skipped
#     downstream -- those were condition dummies. New files do NOT need any
#     skipping: use all columns.)
#   * Sample IDs in the output exactly match the QTL model's sample IDs,
#     because both read the same wide expression TSV + metadata.
#
# INPUT expression matrix: wide TSV, first column `gene`, remaining columns are
# sample IDs (same shape as VST.sampleFiltered.metaConverted.txt). Genes are
# rows, samples are columns; we transpose to samples x genes for PEER/PCA.
#
# OUTPUT layout:
#   <out>/covariates/peer/<CT>/peer_K05.tsv
#   <out>/covariates/expression_pcs/<CT>/exprPC_K05.tsv
#   <out>/covariates/logs/generate_<timestamp>.log
#   <out>/covariates/manifest.tsv
# Each covariate file: first column `sample`, then PEER1..PEERK (or exprPC1..K).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

## --------------------------- tiny arg parser ------------------------------ ##
# Avoids an optparse dependency (not always installed on the cluster).
parse_args <- function(args) {
  opt <- list(
    expression        = NULL,
    metadata          = NULL,
    `cell-type-column`= "celltype",
    `sample-id-column`= "sample",
    `condition-column`= "condition",
    `cell-types`      = NULL,
    `all-cell-types`  = FALSE,
    `k-grid`          = "0,5,10,15",
    `n-hvg`           = "3000",
    cpm               = NULL,           # optional CPM matrix for expressed filter
    `min-cpm`         = "1",            # gene expressed if CPM >= min-cpm ...
    `min-samples`     = "5",            # ... in >= min-samples samples of the CT
    method            = "both",         # peer | expression_pc | both
    `output-dir`      = ".",
    manifest          = "manifest.tsv",
    seed              = "42",
    `peer-lib`        = "",             # optional .libPaths() entry for `peer`
    `peer-condition-cov` = TRUE,        # pass condition dummies to PEER
    `pc-remove-condition` = FALSE,      # regress condition out before expr PCA
    standardize       = TRUE,           # per-gene center+scale before PEER/PCA
    overwrite         = FALSE,
    `dry-run`         = FALSE
  )
  flags_boolean <- c("all-cell-types","overwrite","dry-run","pc-remove-condition")
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (!startsWith(a, "--")) stop("Unexpected argument: ", a)
    key <- sub("^--", "", a)
    # handle --no-<flag> negation for the tri-state defaults-on flags
    if (startsWith(key, "no-")) {
      base <- sub("^no-", "", key)
      if (base %in% names(opt)) { opt[[base]] <- FALSE; i <- i + 1; next }
    }
    if (key %in% flags_boolean) { opt[[key]] <- TRUE; i <- i + 1; next }
    if (key %in% c("peer-condition-cov","standardize")) {
      # allow both bare (--standardize) and valued (--standardize TRUE)
      if (i + 1 <= length(args) && !startsWith(args[[i+1]], "--")) {
        opt[[key]] <- as.logical(args[[i+1]]); i <- i + 2
      } else { opt[[key]] <- TRUE; i <- i + 1 }
      next
    }
    if (!(key %in% names(opt))) stop("Unknown option: --", key)
    if (i + 1 > length(args)) stop("Missing value for --", key)
    opt[[key]] <- args[[i+1]]
    i <- i + 2
  }
  opt
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
  cat("Usage: Rscript generate_covariates.R --expression FILE --metadata FILE \\\n",
      "         [--all-cell-types | --cell-types KRT,MEL,FRB] \\\n",
      "         [--k-grid 0,5,10,15] [--n-hvg 3000] \\\n",
      "         [--method peer|expression_pc|both] [--output-dir DIR] \\\n",
      "         [--seed 42] [--peer-lib PATH] [--no-peer-condition-cov] \\\n",
      "         [--no-standardize] [--overwrite] [--dry-run]\n", sep = "")
  quit(status = 0)
}
opt <- parse_args(args)

log_lines <- character(0)
say <- function(...) {
  msg <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
  cat(msg, "\n", sep = "")
  log_lines[[length(log_lines) + 1]] <<- msg
}

## ------------------------------- validate --------------------------------- ##
if (is.null(opt$expression)) stop("--expression is required")
if (is.null(opt$metadata))   stop("--metadata is required")
stopifnot(file.exists(opt$expression), file.exists(opt$metadata))

method <- match.arg(opt$method, c("both","peer","expression_pc"))
seed   <- as.integer(opt$seed)
n_hvg  <- as.integer(opt$`n-hvg`)
min_cpm     <- as.numeric(opt$`min-cpm`)
min_samples <- as.integer(opt$`min-samples`)
k_grid <- sort(unique(as.integer(strsplit(opt$`k-grid`, ",")[[1]])))
if (any(is.na(k_grid))) stop("--k-grid must be comma-separated integers")
if (any(k_grid < 0))    stop("--k-grid values must be >= 0")

ct_col   <- opt$`cell-type-column`
id_col   <- opt$`sample-id-column`
cond_col <- opt$`condition-column`
out_root <- file.path(opt$`output-dir`, "covariates")

say("==== generate_covariates.R start ====")
say("expression = ", opt$expression)
say("metadata   = ", opt$metadata)
say("method     = ", method, " | k-grid = ", paste(k_grid, collapse=","),
    " | n-hvg = ", n_hvg, " | seed = ", seed)
say("expressed-filter = ", if (is.null(opt$cpm)) "OFF (top-HVG on all genes)"
    else paste0("ON via CPM (>=", min_cpm, " in >=", min_samples, " samples)"))
say("standardize = ", opt$standardize, " | peer-condition-cov = ",
    opt$`peer-condition-cov`, " | pc-remove-condition = ",
    opt$`pc-remove-condition`, " | overwrite = ", opt$overwrite,
    " | dry-run = ", opt$`dry-run`)

## ------------------------------ load inputs ------------------------------- ##
say("Loading metadata ...")
meta <- fread(opt$metadata)
for (col in c(id_col, ct_col, cond_col)) {
  if (!(col %in% names(meta))) stop("Metadata is missing column: ", col)
}
meta <- meta[!is.na(get(id_col)) & !is.na(get(ct_col))]
if (anyDuplicated(meta[[id_col]])) {
  dups <- meta[[id_col]][duplicated(meta[[id_col]])]
  stop("Duplicate sample IDs in metadata: ", paste(head(unique(dups)), collapse=", "))
}

say("Loading expression (header only first, to map samples) ...")
expr <- fread(opt$expression)
gene_col <- names(expr)[1]
say("Expression gene column = '", gene_col, "'; ", nrow(expr), " genes x ",
    ncol(expr) - 1, " sample columns")
expr_samples <- setdiff(names(expr), gene_col)

## optional CPM matrix for the "expressed in this cell type" prefilter
cpm <- NULL; cpm_gene_col <- NULL; cpm_samples <- NULL
if (!is.null(opt$cpm)) {
  stopifnot(file.exists(opt$cpm))
  say("Loading CPM for expressed-gene prefilter: ", opt$cpm)
  cpm <- fread(opt$cpm)
  cpm_gene_col <- names(cpm)[1]
  cpm_samples  <- setdiff(names(cpm), cpm_gene_col)
  say("  CPM: ", nrow(cpm), " genes x ", length(cpm_samples), " samples; ",
      "expressed = CPM >= ", min_cpm, " in >= ", min_samples, " samples of the cell type")
}

## which cell types?
all_cts <- sort(unique(as.character(meta[[ct_col]])))
if (isTRUE(opt$`all-cell-types`)) {
  cts <- all_cts
} else if (!is.null(opt$`cell-types`)) {
  cts <- strsplit(opt$`cell-types`, ",")[[1]]
  bad <- setdiff(cts, all_cts)
  if (length(bad)) stop("Requested cell types not in metadata: ", paste(bad, collapse=", "))
} else {
  stop("Provide --all-cell-types or --cell-types KRT,MEL,FRB")
}
say("Cell types to process: ", paste(cts, collapse = ", "))

## lazy-load peer only if needed
peer_ok <- TRUE
if (method %in% c("both","peer")) {
  if (nzchar(opt$`peer-lib`)) .libPaths(c(opt$`peer-lib`, .libPaths()))
  peer_ok <- requireNamespace("peer", quietly = TRUE)
  if (!peer_ok) {
    say("WARNING: `peer` package not available",
        if (nzchar(opt$`peer-lib`)) paste0(" (checked --peer-lib ", opt$`peer-lib`, ")") else "",
        ". PEER outputs will be SKIPPED; expression_pc still runs if requested.")
  } else {
    suppressMessages(requireNamespace("peer"))
  }
}

## ---------------------------- helper functions ---------------------------- ##
zpad <- function(k) sprintf("K%02d", k)

# per-gene standardization (rows = genes). Returns matrix same orientation.
standardize_genes <- function(m) {
  m <- t(scale(t(m), center = TRUE, scale = TRUE))  # scale each gene
  m
}

# select samples for a cell type; deterministic ordering; loud about drops.
select_samples <- function(ct) {
  meta_ct <- meta[as.character(get(ct_col)) == ct]
  wanted  <- as.character(meta_ct[[id_col]])
  present <- intersect(wanted, expr_samples)
  dropped <- setdiff(wanted, expr_samples)
  extra   <- setdiff(intersect(expr_samples, wanted), wanted)  # always empty; sanity
  if (length(dropped)) {
    say("  NOTE ", ct, ": ", length(dropped), " metadata sample(s) have no ",
        "expression column and are dropped: ",
        paste(head(dropped, 6), collapse=", "),
        if (length(dropped) > 6) " ..." else "")
  }
  # deterministic order: sort sample IDs
  present <- sort(present)
  if (anyDuplicated(present)) stop("Duplicate samples after selection for ", ct)
  list(samples = present, n_dropped = length(dropped))
}

# genes expressed in this cell type's pooled samples (CPM >= min_cpm in
# >= min_samples samples). Mirrors the old get_expressed_genes(). Uses the same
# cell-type sample set as the expression matrix (matched via metadata IDs).
expressed_genes_for_ct <- function(samples) {
  if (is.null(cpm)) return(NULL)
  cpm_use <- intersect(samples, cpm_samples)
  if (length(cpm_use) == 0) {
    say("  WARNING: no CPM columns match this cell type's samples; expressed filter skipped")
    return(NULL)
  }
  if (length(cpm_use) < length(samples)) {
    say("  NOTE: ", length(samples) - length(cpm_use),
        " sample(s) lack a CPM column; expressed filter uses the ",
        length(cpm_use), " with CPM.")
  }
  M <- as.matrix(cpm[, cpm_use, with = FALSE])
  storage.mode(M) <- "double"
  n_expr <- rowSums(M >= min_cpm, na.rm = TRUE)
  cpm[[cpm_gene_col]][n_expr >= min_samples]
}

# build standardized samples x genes matrix for a cell type.
# expressed: optional vector restricting genes BEFORE top-HVG selection.
build_expr_matrix <- function(ct, samples, expressed = NULL) {
  sub <- expr[, c(gene_col, samples), with = FALSE]
  M <- as.matrix(sub[, samples, with = FALSE])
  rownames(M) <- sub[[gene_col]]            # genes x samples
  storage.mode(M) <- "double"
  n_expressed <- NA_integer_
  if (!is.null(expressed)) {
    keep_expr <- rownames(M) %in% expressed
    M <- M[keep_expr, , drop = FALSE]
    n_expressed <- nrow(M)
    say("  expressed genes in ", ct, " (CPM>=", min_cpm, " in >=", min_samples,
        " samples), present in expression matrix: ", n_expressed)
    if (n_expressed < n_hvg)
      say("  NOTE: only ", n_expressed, " expressed genes (< n-hvg ", n_hvg,
          "); all expressed genes used.")
  }
  # drop genes with NA or zero variance
  v <- apply(M, 1, function(x) stats::var(x, na.rm = TRUE))
  keep <- is.finite(v) & v > 0
  M <- M[keep, , drop = FALSE]; v <- v[keep]
  # top n_hvg by variance (among expressed genes, if filtered above)
  ord <- order(v, decreasing = TRUE)
  M <- M[ord[seq_len(min(n_hvg, nrow(M)))], , drop = FALSE]
  n_genes <- nrow(M)
  if (isTRUE(opt$standardize)) M <- standardize_genes(M)
  # a gene can go constant only if input constant (already removed); guard anyway
  M <- M[is.finite(rowSums(M)), , drop = FALSE]
  list(mat = M, n_genes = n_genes, n_expressed = n_expressed)  # genes x samples, standardized
}

write_cov <- function(df_mat, samples, prefix, ct, k, method_dir) {
  # df_mat: samples x K (may be NULL for K=0)
  dir_out <- file.path(out_root, method_dir, ct)
  if (!isTRUE(opt$`dry-run`)) dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
  fpath <- file.path(dir_out, paste0(prefix, "_", zpad(k), ".tsv"))
  if (file.exists(fpath) && !isTRUE(opt$overwrite)) {
    say("  SKIP (exists): ", fpath)
    return(list(path = fpath, written = FALSE, ncol = NA_integer_))
  }
  if (k == 0 || is.null(df_mat)) {
    out <- data.table(sample = samples)
  } else {
    stopifnot(nrow(df_mat) == length(samples))
    out <- data.table(sample = samples)
    cn <- colnames(df_mat)
    for (j in seq_len(ncol(df_mat))) out[[cn[j]]] <- df_mat[, j]
  }
  if (!isTRUE(opt$`dry-run`)) fwrite(out, fpath, sep = "\t")
  say("  wrote ", fpath, "  (", nrow(out), " samples x ", ncol(out) - 1, " covariates)")
  list(path = fpath, written = TRUE, ncol = ncol(out) - 1)
}

## ----------------------------- PEER per-K --------------------------------- ##
fit_peer_K <- function(M_std, samples, ct, k, cond_vec) {
  # M_std: genes x samples (standardized). PEER wants samples x genes.
  Y <- t(M_std)
  model <- peer::PEER()
  peer::PEER_setPhenoMean(model, Y)
  ncov <- 0L
  if (isTRUE(opt$`peer-condition-cov`)) {
    cov <- stats::model.matrix(~ 0 + cond_vec)
    peer::PEER_setCovariates(model, cov)
    ncov <- ncol(cov)
  }
  peer::PEER_setNk(model, k)
  peer::PEER_update(model)
  X <- peer::PEER_getX(model)                 # samples x (ncov + k)  [cov first]
  total <- ncol(X)
  if (total == ncov + k) {
    factors <- X[, (ncov + 1):(ncov + k), drop = FALSE]
  } else if (total == k) {
    factors <- X                              # no covariates prepended
  } else if (total > k) {
    # be defensive: assume learned factors are the trailing k columns
    say("  NOTE PEER getX returned ", total, " cols (expected ", ncov + k,
        "); taking trailing ", k, " as factors.")
    factors <- X[, (total - k + 1):total, drop = FALSE]
  } else {
    stop("PEER returned fewer columns (", total, ") than requested K=", k)
  }
  colnames(factors) <- paste0("PEER", seq_len(k))
  rownames(factors) <- samples
  factors
}

## ------------------------- expression PCs (once) -------------------------- ##
compute_expr_pcs <- function(M_std, samples, maxK, cond_vec = NULL) {
  # M_std: genes x samples (standardized). PCA over samples.
  X <- t(M_std)                               # samples x genes
  if (isTRUE(opt$`pc-remove-condition`) && !is.null(cond_vec) &&
      length(unique(cond_vec)) > 1) {
    # regress condition out of every gene, then PCA on residuals -- so PCs
    # capture hidden/technical structure only, mirroring PEER's condition cov.
    X <- stats::residuals(stats::lm(X ~ cond_vec))
  }
  pc <- stats::prcomp(X, center = TRUE, scale. = FALSE)
  scores <- pc$x[, seq_len(min(maxK, ncol(pc$x))), drop = FALSE]
  # deterministic sign: make the largest-|loading| entry of each PC positive
  for (j in seq_len(ncol(scores))) {
    load_j <- pc$rotation[, j]
    if (load_j[which.max(abs(load_j))] < 0) scores[, j] <- -scores[, j]
  }
  colnames(scores) <- paste0("exprPC", seq_len(ncol(scores)))
  rownames(scores) <- samples
  scores
}

## ------------------------------- main loop -------------------------------- ##
set.seed(seed)
manifest <- data.table()
ts_start <- format(Sys.time(), "%Y%m%d_%H%M%S")

for (ct in cts) {
  say("---- cell type: ", ct, " ----")
  sel <- select_samples(ct)
  samples <- sel$samples
  if (length(samples) < 5) { say("  SKIP ", ct, ": <5 samples"); next }

  expressed <- expressed_genes_for_ct(samples)
  bm <- build_expr_matrix(ct, samples, expressed = expressed)
  M_std <- bm$mat
  n_genes <- bm$n_genes
  say("  ", length(samples), " samples, ", n_genes, " genes used (top ", n_hvg, " HVG)")

  # condition vector aligned to `samples`
  cond_vec <- factor(as.character(meta[[cond_col]][match(samples, meta[[id_col]])]))

  add_manifest <- function(mname, mdir, prefix, k, res) {
    manifest <<- rbind(manifest, data.table(
      cell_type = ct, method = mname, K = k,
      n_samples = length(samples), n_genes = n_genes,
      n_expressed = bm$n_expressed,
      n_covariates = if (is.na(res$ncol)) k else res$ncol,
      output_file = res$path, written = res$written,
      standardized = isTRUE(opt$standardize),
      peer_condition_cov = if (mname == "peer") isTRUE(opt$`peer-condition-cov`) else NA,
      pc_condition_removed = if (startsWith(mname, "expression_pc")) isTRUE(opt$`pc-remove-condition`) else NA,
      input_expression = normalizePath(opt$expression, mustWork = FALSE),
      seed = seed, timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ), fill = TRUE)
    say("  [", mname, " K=", k, "] samples=", length(samples), " genes=", n_genes,
        " -> ", basename(res$path))
  }

  ## ---- PEER ----
  if (method %in% c("both","peer")) {
    for (k in k_grid) {
      if (k == 0) {
        res <- write_cov(NULL, samples, "peer", ct, 0, "peer")
        add_manifest("peer", "peer", "peer", 0, res)
      } else if (!peer_ok) {
        say("  SKIP peer K=", k, " (peer not installed)")
      } else if (isTRUE(opt$`dry-run`)) {
        say("  DRY-RUN would fit PEER K=", k, " for ", ct)
        res <- write_cov(matrix(0, length(samples), k,
                                dimnames=list(NULL, paste0("PEER", seq_len(k)))),
                         samples, "peer", ct, k, "peer")
        add_manifest("peer", "peer", "peer", k, res)
      } else {
        say("  fitting PEER with Nk=", k, " (independent fit) ...")
        factors <- fit_peer_K(M_std, samples, ct, k, cond_vec)
        res <- write_cov(factors, samples, "peer", ct, k, "peer")
        add_manifest("peer", "peer", "peer", k, res)
      }
    }
  }

  ## ---- expression PCs ----
  if (method %in% c("both","expression_pc")) {
    # keep condition-removed PCs in a distinct method/folder so they can coexist
    # with the raw-PC set in the same --output-dir (no filename collision).
    pc_removed  <- isTRUE(opt$`pc-remove-condition`)
    pc_method   <- if (pc_removed) "expression_pc_condremoved" else "expression_pc"
    pc_dir      <- if (pc_removed) "expression_pcs_condremoved" else "expression_pcs"
    maxK <- max(k_grid)
    pcs <- NULL
    if (maxK > 0 && !isTRUE(opt$`dry-run`)) {
      pcs <- compute_expr_pcs(M_std, samples, maxK, cond_vec)   # computed ONCE
    }
    for (k in k_grid) {
      if (k == 0) {
        res <- write_cov(NULL, samples, "exprPC", ct, 0, pc_dir)
      } else if (isTRUE(opt$`dry-run`)) {
        say("  DRY-RUN would compute exprPC K=", k, " for ", ct,
            if (pc_removed) " (condition removed)" else "")
        res <- write_cov(matrix(0, length(samples), k,
                                dimnames=list(NULL, paste0("exprPC", seq_len(k)))),
                         samples, "exprPC", ct, k, pc_dir)
      } else {
        sub <- pcs[, seq_len(min(k, ncol(pcs))), drop = FALSE]   # SUBSET is valid for PCs
        res <- write_cov(sub, samples, "exprPC", ct, k, pc_dir)
      }
      add_manifest(pc_method, pc_dir, "exprPC", k, res)
    }
  }
}

## ------------------------------- outputs ---------------------------------- ##
if (!isTRUE(opt$`dry-run`)) {
  dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
  man_path <- file.path(out_root, opt$manifest)
  # append if exists (idempotent reruns), else create
  if (file.exists(man_path) && !isTRUE(opt$overwrite)) {
    old <- fread(man_path)
    manifest <- unique(rbind(old, manifest, fill = TRUE),
                       by = c("cell_type","method","K"), fromLast = TRUE)
  }
  fwrite(manifest, man_path, sep = "\t")
  say("Manifest written: ", man_path, " (", nrow(manifest), " rows)")

  log_dir <- file.path(out_root, "logs")
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  log_path <- file.path(log_dir, paste0("generate_", ts_start, ".log"))
  writeLines(log_lines, log_path)
  say("Log written: ", log_path)
} else {
  say("DRY-RUN complete: no files written.")
}
say("==== generate_covariates.R done ====")

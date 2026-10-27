# =============================================================================
# load_covariate.R  --  helper to pull one covariate matrix by
# (method, cell_type, K) from a generate_covariates.R output directory.
#
# source() this in the QTL pipeline / optimizer, then:
#   cov <- load_covariate("/pi/.../peer_opt", "peer", "KRT", 10)
#   cov_cols <- setdiff(colnames(cov), "sample")   # PEER1..PEER10 (all of them)
#
# Returns a data.table: first column `sample` (exact QTL sample IDs), then the K
# covariate columns (PEER1..PEERK or exprPC1..exprPCK). K = 0 returns a
# sample-only table (no covariate columns) -- the "no expression covariate"
# baseline. Method is one of:
#   "peer"                       -> peer/<CT>/peer_K<kk>.tsv
#   "expression_pc"              -> expression_pcs/<CT>/exprPC_K<kk>.tsv
#   "expression_pc_condremoved"  -> expression_pcs_condremoved/<CT>/exprPC_K<kk>.tsv
# =============================================================================

suppressPackageStartupMessages(library(data.table))

.cov_method_dir <- c(
  peer                      = "peer",
  expression_pc             = "expression_pcs",
  expression_pc_condremoved = "expression_pcs_condremoved"
)
.cov_prefix <- c(
  peer                      = "peer",
  expression_pc             = "exprPC",
  expression_pc_condremoved = "exprPC"
)

load_covariate <- function(output_dir, method, cell_type, K) {
  if (!(method %in% names(.cov_method_dir)))
    stop("Unknown method '", method, "'. Use one of: ",
         paste(names(.cov_method_dir), collapse = ", "))
  mdir   <- .cov_method_dir[[method]]
  prefix <- .cov_prefix[[method]]
  fname  <- sprintf("%s_K%02d.tsv", prefix, as.integer(K))
  fpath  <- file.path(output_dir, "covariates", mdir, cell_type, fname)
  if (!file.exists(fpath)) stop("Covariate file not found: ", fpath)

  cov <- data.table::fread(fpath)
  data.table::setnames(cov, 1, "sample")           # first col is the sample ID

  n_cov <- ncol(cov) - 1L
  if (n_cov != as.integer(K))
    warning("Expected ", K, " covariate columns in ", fpath, " but found ", n_cov)
  # NOTE: use ALL non-`sample` columns. Do NOT drop PEER1..PEER4 -- these files
  # already contain only learned factors; condition dummies were stripped at
  # generation time.
  cov[]
}

# Convenience: the covariate term string for a model formula, or "" for K = 0.
covariate_terms <- function(cov) {
  cc <- setdiff(colnames(cov), "sample")
  if (length(cc) == 0) return("")
  paste("+", paste(cc, collapse = " + "))
}

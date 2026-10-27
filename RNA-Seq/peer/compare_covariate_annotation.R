#!/usr/bin/env Rscript
# =============================================================================
# compare_covariate_annotation.R
#
# Show how each factor maps onto KNOWN metadata variables ACROSS a grid of K,
# for a given cell type and method. Produces one faceted heatmap per
# (cell_type, method): rows = factor index (1..K), columns = metadata variable,
# fill = R^2 (variance of that factor explained by the variable), facet = K.
#
# This makes it easy to see, e.g., "at K=8 no factor tracks condition, but at
# K=10 two factors do."  NOTE: PEER fits are INDEPENDENT per K (not nested), so
# factor #k at K=10 is not the same vector as factor #k at K=8 -- the facets are
# separate fits. Expression PCs ARE nested (exprPC1..8 of K=8 == first 8 of K=15).
#
# Usage:
#   Rscript compare_covariate_annotation.R --output-dir peer_opt \
#           --metadata metadata.sampleFiltered.txt
#   # optional: --cell-types KRT,MEL  --methods peer  --k-grid 5,8,10,15
#   #           --vars condition,sequencing_platform,phase,total_reads
#   #           --max-levels 60
# =============================================================================

suppressPackageStartupMessages(library(data.table))
if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("ggplot2 is required for this script")
library(ggplot2)

parse_args <- function(args) {
  opt <- list(`output-dir` = NULL, manifest = "manifest.tsv", metadata = NULL,
              `cell-types` = NULL, methods = NULL, `k-grid` = NULL,
              vars = NULL, `sample-id-column` = "sample", `max-levels` = "60")
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!(key %in% names(opt))) stop("Unknown option: --", key)
    opt[[key]] <- args[[i + 1]]; i <- i + 2
  }
  opt
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opt$`output-dir`)) stop("--output-dir is required")
if (is.null(opt$metadata))     stop("--metadata is required")
id_col     <- opt$`sample-id-column`
max_levels <- as.integer(opt$`max-levels`)

meta <- fread(opt$metadata)
man_path <- file.path(opt$`output-dir`, "covariates", opt$manifest)
if (!file.exists(man_path)) man_path <- file.path(opt$`output-dir`, opt$manifest)
if (!file.exists(man_path)) stop("Manifest not found under --output-dir")
man <- fread(man_path)[K > 0]

if (!is.null(opt$`cell-types`)) man <- man[cell_type %in% strsplit(opt$`cell-types`, ",")[[1]]]
if (!is.null(opt$methods))      man <- man[method    %in% strsplit(opt$methods, ",")[[1]]]
if (!is.null(opt$`k-grid`))     man <- man[K %in% as.integer(strsplit(opt$`k-grid`, ",")[[1]])]
if (nrow(man) == 0) stop("No matching covariate files in manifest after filtering")

cand_vars <- setdiff(names(meta), c(id_col, "celltype"))
if (!is.null(opt$vars)) cand_vars <- strsplit(opt$vars, ",")[[1]]

# Resolve each covariate file relative to --output-dir when the manifest's
# absolute path (recorded on the cluster) doesn't exist locally -- so the folder
# can be moved/copied anywhere and this still works.
method_subdir <- c(peer = "peer", expression_pc = "expression_pcs",
                   expression_pc_condremoved = "expression_pcs_condremoved")
resolve_path <- function(output_file, method, cell_type) {
  if (file.exists(output_file)) return(output_file)
  alt <- file.path(opt$`output-dir`, "covariates",
                   method_subdir[[method]], cell_type, basename(output_file))
  if (file.exists(alt)) return(alt)
  output_file  # return original; caller reports MISSING
}
man[, resolved := mapply(resolve_path, output_file, method, cell_type)]

# association of one factor vector against one metadata variable -> R^2
assoc_r2 <- function(y, xv, is_num) {
  dat <- data.frame(y = y, x = if (is_num) as.numeric(xv) else factor(xv))
  dat <- dat[stats::complete.cases(dat), ]
  if (nrow(dat) < 3 || length(unique(dat$x)) < 2) return(c(r2 = NA, p = NA))
  fit <- tryCatch(stats::lm(y ~ x, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(c(r2 = NA, p = NA))
  s <- summary(fit); fst <- s$fstatistic
  p <- if (!is.null(fst)) stats::pf(fst[1], fst[2], fst[3], lower.tail = FALSE) else NA
  c(r2 = s$r.squared, p = as.numeric(p))
}

annotate_file <- function(f, K) {
  if (!file.exists(f)) { message("  MISSING: ", f); return(NULL) }
  cov <- fread(f)
  fac_cols <- setdiff(names(cov), "sample")
  m <- meta[match(cov$sample, get(id_col))]
  keep <- !is.na(m[[id_col]]); cov <- cov[keep]; m <- m[keep]
  rows <- list()
  for (v in cand_vars) {
    if (!(v %in% names(m))) next
    x <- m[[v]]; if (all(is.na(x))) next
    is_num <- is.numeric(x); nlev <- length(unique(x[!is.na(x)]))
    if (nlev < 2) next
    if (!is_num && max_levels > 0 && nlev > max_levels) next
    for (fi in seq_along(fac_cols)) {
      rr <- assoc_r2(as.numeric(cov[[fac_cols[fi]]]), x, is_num)
      rows[[length(rows) + 1]] <- data.table(K = K, factor_idx = fi,
        factor = fac_cols[fi], variable = v, r2 = rr[["r2"]], p = rr[["p"]])
    }
  }
  rbindlist(rows)
}

out_dir <- file.path(opt$`output-dir`, "covariates")
combos <- unique(man[, .(cell_type, method)])
for (ci in seq_len(nrow(combos))) {
  ct <- combos$cell_type[ci]; mth <- combos$method[ci]
  sub <- man[cell_type == ct & method == mth][order(K)]
  message("== ", ct, " / ", mth, "  (K = ", paste(sub$K, collapse=","), ") ==")
  dat <- rbindlist(lapply(seq_len(nrow(sub)),
                          function(i) annotate_file(sub$resolved[i], sub$K[i])))
  if (is.null(dat) || nrow(dat) == 0) { message("  nothing to plot"); next }

  dat[, Klab := factor(paste0("K=", K), levels = paste0("K=", sort(unique(K))))]
  dat[, sig := ifelse(is.na(p), "", ifelse(p < 1e-4, "***",
               ifelse(p < 1e-3, "**", ifelse(p < 0.05, "*", ""))))]
  maxK <- max(dat$factor_idx)
  prefix <- if (mth == "peer") "PEER" else "exprPC"

  g <- ggplot(dat, aes(x = variable, y = factor_idx, fill = r2)) +
    geom_tile(color = "grey88") +
    geom_text(aes(label = sig), size = 2.6) +
    facet_wrap(~ Klab, nrow = 1) +
    scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, 1),
                        name = expression(R^2)) +
    scale_y_reverse(breaks = seq_len(maxK),
                    labels = paste0(prefix, seq_len(maxK))) +
    labs(title = paste0(ct, " - ", mth, ": factor vs. known metadata across K"),
         subtitle = "R^2 of each factor explained by metadata (* p<.05, ** p<1e-3, *** p<1e-4)",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())

  n_vars <- length(unique(dat$variable)); n_K <- length(unique(dat$K))
  out_pdf <- file.path(out_dir, sprintf("annotation_compare_%s_%s.pdf", ct, mth))
  ggsave(out_pdf, g, width = max(7, (0.32 * n_vars + 1.2) * n_K),
         height = max(4, 0.28 * maxK + 2), limitsize = FALSE)
  fwrite(dat[order(K, factor_idx, -r2)],
         file.path(out_dir, sprintf("annotation_compare_%s_%s.tsv", ct, mth)), sep = "\t")
  message("  wrote ", out_pdf)
}
message("done.")

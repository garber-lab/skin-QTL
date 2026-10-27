#!/usr/bin/env Rscript
# =============================================================================
# plot_scree.R  --  scree plot from a plink (--pca) .eigenval file.
#
# Reads the eigenvalues, estimates the noise floor (mean of the smallest few
# PCs), and classifies each PC as real / marginal / noise by its ratio to the
# floor. Writes a PDF scree and prints a table + a suggested number of PCs to
# keep as genotype covariates (smallest set clearly above the floor).
#
# The floor + ratio-to-floor framing is used (rather than % variance) because a
# `--pca K` run only reports the top K eigenvalues, so the true total variance
# is unknown -- the excess above the noise floor is the robust signal.
#
# Usage:
#   Rscript plot_scree.R --eigenval imputed_pca.eigenval [--out scree.pdf]
#           [--floor-n 5] [--real 1.15] [--marginal 1.05] [--label PCA]
# =============================================================================
suppressPackageStartupMessages(library(data.table))

parse_args <- function(a) {
  o <- list(eigenval=NULL, out=NULL, `floor-n`="5", real="1.25", marginal="1.10", label="genotype PCA")
  i <- 1; while (i <= length(a)) { k <- sub("^--","",a[[i]])
    if (!(k %in% names(o))) stop("Unknown option --", k); o[[k]] <- a[[i+1]]; i <- i+2 }; o
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opt$eigenval)) stop("--eigenval is required")
if (is.null(opt$out)) opt$out <- sub("\\.eigenval$", "_scree.pdf", opt$eigenval)
floor_n  <- as.integer(opt$`floor-n`)
thr_real <- as.numeric(opt$real)
thr_marg <- as.numeric(opt$marginal)

ev <- fread(opt$eigenval, header = FALSE)[[1]]
ev <- ev[is.finite(ev)]
K <- length(ev)
if (K < 3) stop("Need >=3 eigenvalues")
floor_n <- min(floor_n, K - 1)
floor <- mean(tail(sort(ev, decreasing = TRUE), floor_n))   # mean of the smallest floor_n
ratio <- ev / floor

cls <- ifelse(ratio >= thr_real, "real", ifelse(ratio >= thr_marg, "marginal", "noise"))
n_real <- sum(cls == "real"); n_marg <- sum(cls == "marginal")

tab <- data.table(PC = paste0("PC", seq_len(K)), eigenvalue = round(ev, 4),
                  ratio_to_floor = round(ratio, 3), class = cls)
cat("Noise floor (mean of smallest ", floor_n, " PCs): ", round(floor, 4), "\n", sep = "")
print(tab)
cat("\nSuggested genotype PCs to keep: ", n_real,
    " clearly-real", if (n_marg) paste0(" (+", n_marg, " marginal -> test up to ", n_real + n_marg, ")") else "",
    "\nPick the smallest nGPC where QTL lambda ~ 1.\n", sep = "")

## ---- scree PDF ----
col_map <- c(real = "#185FA5", marginal = "#EF9F27", noise = "#888780")
pdf(opt$out, width = 7, height = 5)
plot(seq_len(K), ev, type = "b", pch = 16, cex = 1.3, col = col_map[cls],
     xaxt = "n", xlab = NA, ylab = "eigenvalue",
     ylim = c(min(ev) * 0.95, max(ev) * 1.03),
     main = paste0("Scree: ", opt$label))
axis(1, at = seq_len(K), labels = tab$PC, las = 2)
abline(h = floor, lty = 2, col = "#993C1D")
text(K, floor, "noise floor", pos = 3, col = "#993C1D", cex = 0.8)
if (n_real > 0) abline(v = n_real + 0.5, lty = 3, col = "grey60")
legend("topright", c("real", "marginal", "noise", "floor"),
       col = c(col_map, "#993C1D"), pch = c(16,16,16,NA), lty = c(NA,NA,NA,2), bty = "n", cex = 0.9)
dev.off()
cat("Wrote scree plot: ", opt$out, "\n", sep = "")

#!/usr/bin/env Rscript
# =============================================================================
# format_eigenvec_for_qtl.R
#
# Convert a plink2 --pca .eigenvec into the EXACT layout the QTL pipeline and
# the covariate optimizer expect (matching genotype_PCs_for_modeling_07242025.txt):
#
#   FID              donor   donor_num   PC1        PC2 ...
#   skineQTL-F117    F117    117         -0.048...  ...
#
# step1_run_QTL_chunk.R / opt_fit_qtl_subset.R merge genotype PCs by `donor`, and
# take pc_cols = setdiff(names, c("FID","donor","donor_num")). So the column names
# and the `donor` values must match your metadata donor IDs (e.g. F117, F75).
#
# donor is derived from the sample IID by stripping the FID prefix (default
# "skineQTL-") and a single leading zero after F (F075 -> F75), matching how the
# QTL genotype tables clean donor IDs. Verify the printed examples line up with
# your metadata before using downstream.
#
# Usage:
#   Rscript format_eigenvec_for_qtl.R --eigenvec imputed_pca.eigenvec \
#           --out genotype_PCs_for_modeling_imputed.txt [--strip-prefix skineQTL-]
# =============================================================================
suppressPackageStartupMessages(library(data.table))

parse_args <- function(a) {
  o <- list(eigenvec=NULL, out=NULL, `strip-prefix`="skineQTL-")
  i <- 1; while (i <= length(a)) { k <- sub("^--","",a[[i]])
    if (!(k %in% names(o))) stop("Unknown option --", k); o[[k]] <- a[[i+1]]; i <- i+2 }; o
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opt$eigenvec) || is.null(opt$out)) stop("--eigenvec and --out are required")

ev <- fread(opt$eigenvec)
setnames(ev, sub("^#", "", names(ev)))          # plink2 header starts with '#FID' or '#IID'

# sample id column: prefer IID, else the first column
iid_col <- if ("IID" %in% names(ev)) "IID" else names(ev)[1]
pc_cols <- grep("^PC[0-9]+$", names(ev), value = TRUE)
if (!length(pc_cols)) stop("No PC columns found in ", opt$eigenvec)

iid <- as.character(ev[[iid_col]])
donor <- sub(paste0("^", opt$`strip-prefix`), "", iid)   # strip FID prefix
donor <- sub("^F0", "F", donor)                          # F075 -> F75 (single leading zero)
donor_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", donor)))

out <- data.table(FID = iid, donor = donor, donor_num = donor_num)
out <- cbind(out, ev[, ..pc_cols])

if (anyDuplicated(out$donor)) {
  dups <- out$donor[duplicated(out$donor)]
  stop("Duplicate donor IDs after cleaning: ", paste(unique(dups), collapse=", "),
       " -- check --strip-prefix / naming.")
}

fwrite(out, opt$out, sep = "\t")
cat("Wrote ", opt$out, ": ", nrow(out), " donors x ", length(pc_cols), " PCs\n", sep = "")
cat("Verify these donor IDs match your metadata:\n")
print(utils::head(out[, .(FID, donor, donor_num)], 5))

#!/usr/bin/env Rscript
# =============================================================================
# assign_ancestry.R
#
# From a combined study+1000G PCA, label each study donor with a superpopulation
# by k-nearest-neighbours on the reference samples (their known SuperPop), using
# the top PCs. Flags likely-admixed donors (no clear majority among neighbours).
#
# Output table: donor, predicted_superpop, confidence (top vote fraction),
#               and per-superpop vote fractions; plus a PC1-PC2 overlay PDF.
#
# Usage:
#   Rscript assign_ancestry.R --eigenvec combined_pca.eigenvec \
#     --ref-labels ref_superpop.txt --study-iids study_iids.txt \
#     --out ancestry_estimation.txt [--plot ancestry_PCA_1kg.pdf] \
#     [--k 10] [--npc 6] [--strip-prefix skineQTL-]
# =============================================================================
suppressPackageStartupMessages(library(data.table))

parse_args <- function(a) {
  o <- list(eigenvec=NULL, `ref-labels`=NULL, `study-iids`=NULL, out=NULL,
            plot=NULL, k="10", npc="6", `strip-prefix`="skineQTL-")
  i <- 1; while (i <= length(a)) { k <- sub("^--","",a[[i]])
    if (!(k %in% names(o))) stop("Unknown option --", k); o[[k]] <- a[[i+1]]; i <- i+2 }; o
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
for (r in c("eigenvec","ref-labels","study-iids","out"))
  if (is.null(opt[[r]])) stop("--", r, " is required")
K   <- as.integer(opt$k)
NPC <- as.integer(opt$npc)

ev <- fread(opt$eigenvec); setnames(ev, sub("^#", "", names(ev)))
iid_col <- if ("IID" %in% names(ev)) "IID" else names(ev)[1]
pc_cols <- grep("^PC[0-9]+$", names(ev), value = TRUE)[seq_len(NPC)]
ev[, iid := as.character(get(iid_col))]

reflab <- fread(opt$`ref-labels`, header = FALSE, col.names = c("iid","superpop"))
study_iids <- unique(trimws(readLines(opt$`study-iids`)))

ev[, is_study := iid %in% study_iids]
ev <- merge(ev, reflab, by = "iid", all.x = TRUE)

ref <- ev[is_study == FALSE & !is.na(superpop)]
stu <- ev[is_study == TRUE]
if (nrow(ref) == 0) stop("No reference samples with labels matched the eigenvec IIDs.")
if (nrow(stu) == 0) stop("No study samples matched --study-iids in the eigenvec.")
cat("Reference:", nrow(ref), "labelled samples |", nrow(stu), "study donors |",
    "SuperPops:", paste(sort(unique(ref$superpop)), collapse=","), "\n")

Rm <- as.matrix(ref[, ..pc_cols]); Sm <- as.matrix(stu[, ..pc_cols])
pops <- sort(unique(ref$superpop))

# k-NN vote for each study sample
res <- rbindlist(lapply(seq_len(nrow(Sm)), function(i) {
  d <- sqrt(colSums((t(Rm) - Sm[i, ])^2))          # euclidean in PC space
  nn <- ref$superpop[order(d)[seq_len(K)]]
  fr <- prop.table(table(factor(nn, levels = pops)))
  top <- names(which.max(fr))
  out <- data.table(iid = stu$iid[i], predicted_superpop = top,
                    confidence = round(max(fr), 3))
  for (p in pops) out[[paste0("vote_", p)]] <- round(as.numeric(fr[p]), 3)
  out
}))

# clean donor id to match metadata, flag admixed (no strong majority)
res[, donor := sub(paste0("^", opt$`strip-prefix`), "", iid)]
res[, donor := sub("^F0", "F", donor)]
res[, admixed_flag := confidence < 0.6]
setcolorder(res, c("donor","iid","predicted_superpop","confidence","admixed_flag"))
fwrite(res, opt$out, sep = "\t")
cat("Wrote", opt$out, "\n")
cat("\nSuperpopulation assignment (n donors):\n"); print(table(res$predicted_superpop))
cat("Likely admixed (confidence < 0.6):", sum(res$admixed_flag), "\n")

## ---- optional overlay plot (PC1 vs PC2) ----
if (!is.null(opt$plot)) {
  cols <- c(AFR="#D85A30", AMR="#639922", EAS="#EF9F27", EUR="#534AB7", SAS="#378ADD")
  pdf(opt$plot, width = 7, height = 6)
  plot(ref$PC1, ref$PC2, col = cols[ref$superpop], pch = 16, cex = 0.5,
       xlab = "PC1", ylab = "PC2", main = "Study donors on 1000G ancestry PCA")
  points(stu$PC1, stu$PC2, pch = 21, bg = "black", col = "white", cex = 1.1)
  legend("topright", c(names(cols), "study"), pch = c(rep(16,5),21),
         col = c(cols, "black"), pt.bg = c(rep(NA,5),"black"), bty = "n", cex = 0.8)
  dev.off()
  cat("Wrote plot", opt$plot, "\n")
}

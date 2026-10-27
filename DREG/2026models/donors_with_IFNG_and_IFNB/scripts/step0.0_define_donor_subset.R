#!/usr/bin/env Rscript
# =============================================================================
# step0.0_define_donor_subset.R  --  define the matched-donor subset.
#
# Subset = donors that
#   (a) have RNA-seq expression in BOTH IFNG and IFNB in the given celltype, AND
#   (b) have genotype data.
# One file written per celltype:
#   <rootdir>/data/donor_subset_<ct>.txt   (one donor id per line)
#
# Rationale: the genome-wide pipeline fits each condition on whatever donors are
# available, so IFNG and IFNB rest on different (and different-sized) donor sets.
# Restricting to this subset makes IFNG vs IFNB a matched-donor comparison, so an
# IFNB eGene deficit cannot be a donor-count / detection-power artifact.
#
# NOTE: the subset constrains WHICH donors enter every model, not which
# conditions. Those donors' PBS / TNF samples are still used (PBS is needed as the
# reQTL / condition reference); the point is only that IFNG and IFNB share one
# donor set.
#
# Usage:
#   Rscript step0.0_define_donor_subset.R <rootdir> <META_FILE> <GENO_DONORS_FILE>
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
rootdir          <- if (length(args) >= 1) args[1] else "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB"
META_FILE        <- if (length(args) >= 2) args[2] else "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
GENO_DONORS_FILE <- if (length(args) >= 3) args[3] else paste0(rootdir, "/data/genotyped_donors.txt")

## ---- load inputs ----
meta <- fread(META_FILE)
stopifnot(all(c("sample","donor","condition","celltype") %in% names(meta)))

# genotyped donor ids, cleaned to the SAME namespace as meta$donor upstream
# (see step0 shell: bcftools query -l | cut -d'_' -f1 | sed skineQTL- | sed F0->F)
geno_donors <- unique(trimws(readLines(GENO_DONORS_FILE)))
geno_donors <- geno_donors[nzchar(geno_donors)]
message("Genotyped donors: ", length(geno_donors))

## ---- build per-celltype subset ----
celltypes <- c("MEL","KRT","FRB")
summary_rows <- list()

for (ct in celltypes) {
  m <- meta[celltype == ct]

  ifng_donors <- unique(m[condition == "IFNG", donor])
  ifnb_donors <- unique(m[condition == "IFNB", donor])
  both        <- intersect(ifng_donors, ifnb_donors)          # expr in BOTH IFNG & IFNB
  subset      <- sort(intersect(both, geno_donors))           # ... and genotyped

  outfile <- paste0(rootdir, "/data/donor_subset_", ct, ".txt")
  writeLines(as.character(subset), outfile)

  cat(sprintf("[%s] %s: IFNG=%d IFNB=%d both=%d both&genotyped=%d -> %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ct,
              length(ifng_donors), length(ifnb_donors),
              length(both), length(subset), outfile))

  summary_rows[[ct]] <- data.table(
    celltype = ct,
    n_IFNG = length(ifng_donors), n_IFNB = length(ifnb_donors),
    n_both = length(both), n_subset = length(subset)
  )
}

summary_dt <- rbindlist(summary_rows)
fwrite(summary_dt, paste0(rootdir, "/data/donor_subset_summary.txt"), sep = "\t")
message("Wrote per-celltype donor subsets + donor_subset_summary.txt")

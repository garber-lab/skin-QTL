#!/usr/bin/env Rscript
# =============================================================================
# opt_make_folds.R  --  build gene-subset folds for covariate optimization
#
# Partitions the tested-gene universe into N disjoint random folds of ~M genes,
# and writes, per fold, a pairs file + a genotype file in the SAME format as the
# production chunks (pairs_chunk_<tag>.tsv / genotype_pairs_chunk_<tag>.tsv), so
# the fitter and steps 2-4 run per fold with chunk_id = fold tag.
#
# Folds are shared across cell types (the chunks are); each cell type's a-priori
# BH gene family is applied later in step4, so a fold is just a random gene set.
#
# Usage:
#   Rscript opt_make_folds.R --chunks-dir DIR --out-dir DIR \
#           [--n-folds 5] [--n-genes-per-fold 3000] [--seed 42] \
#           [--gene-pool-file genes.txt]   # optional: restrict universe
# Output:
#   <out-dir>/pairs_chunk_fold1.tsv ... genotype_pairs_chunk_fold1.tsv ...
#   <out-dir>/fold_genes_fold1.txt ...
#   <out-dir>/folds_manifest.tsv
# =============================================================================
suppressPackageStartupMessages(library(data.table))

parse_args <- function(a) {
  o <- list(`chunks-dir`=NULL, `out-dir`=NULL, `n-folds`="5",
            `n-genes-per-fold`="3000", seed="42", `gene-pool-file`=NULL,
            `seed-genes-file`=NULL, `n-seed-per-fold`="100",
            `subchunk-genes`="50", vst=NULL, overwrite=FALSE)
  i <- 1; while (i <= length(a)) { k <- sub("^--","",a[[i]])
    if (!(k %in% names(o))) stop("Unknown option --", k)
    if (k == "overwrite") { o[[k]] <- TRUE; i <- i+1; next }   # boolean flag
    o[[k]] <- a[[i+1]]; i <- i+2 }; o
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opt$`chunks-dir`) || is.null(opt$`out-dir`))
  stop("--chunks-dir and --out-dir are required")
n_folds <- as.integer(opt$`n-folds`)
n_per   <- as.integer(opt$`n-genes-per-fold`)
seed    <- as.integer(opt$seed)
n_seed  <- as.integer(opt$`n-seed-per-fold`)
subchunk_genes <- as.integer(opt$`subchunk-genes`)
overwrite <- isTRUE(opt$overwrite)
dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)

# skip-if-exists writers (unless --overwrite); return TRUE if written, FALSE if skipped
fwrite_skip <- function(dt, path) {
  if (file.exists(path) && !overwrite) return(FALSE)
  fwrite(dt, path, sep = "\t"); TRUE
}
writeLines_skip <- function(x, path) {
  if (file.exists(path) && !overwrite) return(FALSE)
  writeLines(x, path); TRUE
}

pair_files <- sort(list.files(opt$`chunks-dir`, pattern = "^pairs_chunk_.*\\.tsv$", full.names = TRUE))
geno_files <- sort(list.files(opt$`chunks-dir`, pattern = "^genotype_pairs_chunk_.*\\.tsv$", full.names = TRUE))
if (!length(pair_files)) stop("No pairs_chunk_*.tsv in ", opt$`chunks-dir`)
if (!length(geno_files)) stop("No genotype_pairs_chunk_*.tsv in ", opt$`chunks-dir`)
message("Found ", length(pair_files), " pair chunks, ", length(geno_files), " genotype chunks")

## ---- load all pairs (keep every column so downstream positional cols work) ----
pairs <- rbindlist(lapply(pair_files, fread), use.names = TRUE, fill = TRUE)
stopifnot("gene_name" %in% names(pairs), "SNP_ID" %in% names(pairs))
message("Total pairs: ", nrow(pairs), " | unique genes: ", uniqueN(pairs$gene_name))

## ---- gene universe ----
pool <- unique(pairs$gene_name)
if (!is.null(opt$`gene-pool-file`)) {
  keep <- unique(trimws(readLines(opt$`gene-pool-file`)))
  pool <- intersect(pool, keep)
  message("Restricted gene pool to ", length(pool), " genes via ", opt$`gene-pool-file`)
}
## ---- optional seed set: guarantee known eGenes in each fold ----
# Seeding a fixed number of previously-significant eGenes per fold gives every
# fold a known-positive set, so the scorer can report RECOVERY (sensitivity)
# per covariate setting, not just raw discovery. lambda is unaffected: it is the
# median chi-sq over millions of cis pairs, and even within a seeded eGene most
# SNPs are null, so ~100 extra eGene genes barely move the median.
seed_genes <- character(0)
if (!is.null(opt$`seed-genes-file`)) {
  sg <- unique(trimws(readLines(opt$`seed-genes-file`)))
  seed_genes <- intersect(sg, pool)                 # only seed testable genes
  message("Seed eGenes: ", length(sg), " provided, ", length(seed_genes),
          " present in gene pool")
}
set.seed(seed)
seed_shuf <- sample(seed_genes)

# disjoint per-fold seed count, capped by availability
n_seed_per <- if (length(seed_genes) == 0) 0L else min(n_seed, floor(length(seed_genes) / n_folds))
if (length(seed_genes) > 0 && n_seed_per < n_seed)
  message("NOTE: only ", length(seed_genes), " seed eGenes for ", n_folds,
          " folds; using ", n_seed_per, "/fold (requested ", n_seed, ").")

# Only the genes actually USED as guaranteed seeds are held out of the fill.
# The fill is a random draw from the REST of the genome -- which still contains
# the non-used eGenes -- so the background carries eGenes at the natural rate
# (representative sample), not an artificial eGene-free complement. The seeded
# 100 give a stable recovery denominator on top of that.
used_seeds <- if (n_seed_per > 0) seed_shuf[seq_len(n_seed_per * n_folds)] else character(0)
fill_pool  <- setdiff(pool, used_seeds)
set.seed(seed + 1L)                    # independent shuffle for the fill
fill_shuf  <- sample(fill_pool)

folds <- vector("list", n_folds); seeds_in_fold <- vector("list", n_folds)
s_cur <- 0L; f_cur <- 0L
for (fi in seq_len(n_folds)) {
  s_i <- if (n_seed_per > 0) seed_shuf[(s_cur + 1L):(s_cur + n_seed_per)] else character(0)
  s_cur <- s_cur + n_seed_per
  n_fill <- n_per - length(s_i)
  fill_i <- fill_shuf[(f_cur + 1L):(f_cur + n_fill)]
  fill_i <- fill_i[!is.na(fill_i)]                  # in case pool runs short
  f_cur  <- f_cur + n_fill
  folds[[fi]] <- c(s_i, fill_i)
  seeds_in_fold[[fi]] <- s_i
}
if (f_cur > length(fill_shuf))
  message("NOTE: fill pool exhausted; some folds are smaller than ", n_per, ".")

## ---- load all genotype rows once, keyed by ID (dedup across chunks) ----
geno <- rbindlist(lapply(geno_files, fread), use.names = TRUE, fill = TRUE)
geno <- geno[!duplicated(ID)]
stopifnot("ID" %in% names(geno))

## ---- optional: load the full VST once, to write tiny per-sub-chunk subsets ----
# so each fit job reads a ~few-gene file instead of the 409MB whole-VST fread
# (drops per-job memory from ~3GB to a few MB and removes the repeated big read).
vst <- NULL; vst_gene_col <- NULL
if (!is.null(opt$vst)) {
  message("Loading VST for per-sub-chunk subsets: ", opt$vst)
  vst <- fread(opt$vst)
  vst_gene_col <- names(vst)[1]                 # first column is the gene id
  message("  VST: ", nrow(vst), " genes x ", ncol(vst) - 1, " samples")
}

manifest <- data.table(); subchunk_index <- data.table()
for (fi in seq_along(folds)) {
  tag <- paste0("fold", fi)
  fg  <- folds[[fi]]
  s_i <- seeds_in_fold[[fi]]

  writeLines_skip(fg, file.path(opt$`out-dir`, paste0("fold_genes_", tag, ".txt")))
  # seed eGene list = the guaranteed recovery set (used as the scorer's denominator)
  writeLines_skip(s_i, file.path(opt$`out-dir`, paste0("fold_seed_egenes_", tag, ".txt")))
  # per-gene annotation: prior_eGene by the FULL known-eGene list (so natural
  # background eGenes in the fill are labelled too); seeded flags the recovery set.
  fwrite_skip(data.table(gene = fg, prior_eGene = fg %in% seed_genes, seeded = fg %in% s_i),
              file.path(opt$`out-dir`, paste0("fold_gene_annot_", tag, ".tsv")))
  n_prior_in_fold <- sum(fg %in% seed_genes)

  # GENE-AWARE sub-chunking: split the fold's genes into groups of
  # `subchunk_genes` (each group keeps every SNP for its genes together, so
  # eigenMT per gene is intact), and write production-sized pair/genotype files.
  # Each sub-chunk fits like a normal chunk (small geno table) and is one array
  # element; results are aggregated per fold before eigenMT->BH.
  grp <- split(fg, ceiling(seq_along(fg) / subchunk_genes))
  fold_pairs <- 0L
  for (si in seq_along(grp)) {
    sub_tag <- sprintf("%s_s%03d", tag, si)
    gsub  <- grp[[si]]
    psub  <- pairs[gene_name %in% gsub]
    snpsub<- unique(psub$SNP_ID)
    gtsub <- geno[ID %in% snpsub]
    fwrite_skip(psub, file.path(opt$`out-dir`, paste0("pairs_chunk_", sub_tag, ".tsv")))
    fwrite_skip(gtsub,file.path(opt$`out-dir`, paste0("genotype_pairs_chunk_", sub_tag, ".tsv")))
    # tiny per-sub-chunk VST subset (rows = this sub-chunk's genes); the fitter
    # reads vst_chunk_<sub_tag>.tsv instead of the whole VST when present.
    if (!is.null(vst)) {
      vsub <- vst[get(vst_gene_col) %in% gsub]
      fwrite_skip(vsub, file.path(opt$`out-dir`, paste0("vst_chunk_", sub_tag, ".tsv")))
    }
    subchunk_index <- rbind(subchunk_index, data.table(
      fold = tag, subchunk = sub_tag, n_genes = length(gsub),
      n_pairs = nrow(psub), n_snps = length(snpsub)))
    fold_pairs <- fold_pairs + nrow(psub)
  }
  # list of sub-chunk tags for this fold (driver reads this to build the fit array)
  writeLines_skip(subchunk_index[fold == tag, subchunk],
                  file.path(opt$`out-dir`, paste0("subchunks_", tag, ".txt")))

  manifest <- rbind(manifest, data.table(
    fold = tag, n_genes = length(fg), n_seed_egenes = length(s_i),
    n_prior_egenes_total = n_prior_in_fold, n_subchunks = length(grp), n_pairs = fold_pairs))
  message(sprintf("  %s: %d genes (%d seeded + %d natural prior eGenes), %d sub-chunks, %d pairs",
                  tag, length(fg), length(s_i), n_prior_in_fold - length(s_i), length(grp), fold_pairs))
}
fwrite(manifest, file.path(opt$`out-dir`, "folds_manifest.tsv"), sep = "\t")
fwrite(subchunk_index, file.path(opt$`out-dir`, "subchunk_index.tsv"), sep = "\t")
message("Wrote ", n_folds, " folds (sub-chunked at ", subchunk_genes,
        " genes) to ", opt$`out-dir`, "  (seed=", seed, ")")

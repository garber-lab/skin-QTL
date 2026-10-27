#!/usr/bin/env Rscript
# =============================================================================
# step0.4_chunk_vst.R  --  chunk the VST expression matrix to match gene chunks.
#
# For each gene_chunk_<id> (already made by step0.2), write vst_chunk_<id>.tsv
# containing only that chunk's gene rows (ALL sample columns kept -- step1 does
# the per-cell-type column select). step1 then reads the tiny per-chunk VST
# instead of the full ~409MB matrix: far less memory per job and no repeated
# whole-VST read. Does NOT touch the gene or genotype chunks.
#
# Reads the full VST ONCE, then writes all subsets. Skips any vst_chunk that
# already exists (idempotent).
#
# Usage:
#   Rscript step0.4_chunk_vst.R <VST_FILE> <CHUNK_DIR>
# =============================================================================
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript step0.4_chunk_vst.R <VST_FILE> <CHUNK_DIR>")
VST_FILE  <- args[1]
CHUNK_DIR <- args[2]

vst <- fread(VST_FILE)
gene_col <- names(vst)[1]                       # first column is the gene id
gfiles <- sort(list.files(CHUNK_DIR, pattern = "^gene_chunk_", full.names = TRUE))
if (!length(gfiles)) stop("No gene_chunk_* files in ", CHUNK_DIR)
cat("VST:", nrow(vst), "genes x", ncol(vst) - 1, "samples |",
    length(gfiles), "gene chunks\n")

n_written <- 0L; n_skipped <- 0L
for (g in gfiles) {
  id  <- sub("^gene_chunk_", "", basename(g))
  out <- file.path(CHUNK_DIR, paste0("vst_chunk_", id, ".tsv"))
  if (file.exists(out)) { n_skipped <- n_skipped + 1L; next }
  genes <- unique(trimws(readLines(g)))
  fwrite(vst[get(gene_col) %in% genes], out, sep = "\t")
  n_written <- n_written + 1L
}
cat("Wrote", n_written, "vst_chunk files; skipped", n_skipped, "existing.\n")

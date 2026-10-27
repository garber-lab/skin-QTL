# decide gene chunks based on SNP:gene pairs
# MATCHED-DONOR variant: parametrized (pairs file, chunk dir, genes/chunk, max pairs)
# so it can be reused for the union-gene set at 25 genes/chunk.
#
# Usage:
#   Rscript step0.2_compile_gene_chunks.R <pairs_file> <chunk_dir> [genes_per_chunk=25] [max_pairs=100000]
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript step0.2_compile_gene_chunks.R <pairs_file> <chunk_dir> [genes_per_chunk=25] [max_pairs=100000]")
}
pairs_file      <- args[[1]]
chunk_dir       <- args[[2]]
genes_per_chunk <- if (length(args) >= 3) as.integer(args[[3]]) else 25L
max_pairs       <- if (length(args) >= 4) as.integer(args[[4]]) else 100000L

# read the (union-filtered) variant-gene pairs
f <- fread(pairs_file)

# summarize number of SNPs per gene
summary <- f[, .(nSNP = .N), by = gene_name]
setnames(summary, "gene_name", "gene")
rm(f)

set.seed(42)
summary <- summary[sample(.N)]

chunks <- list()
current_genes <- character(0)
current_pairs <- 0
chunk_i <- 1

for (k in seq_len(nrow(summary))) {
  g <- summary$gene[k]
  n <- summary$nSNP[k]

  # single gene alone exceeds the cap -- can't split a gene, give it its own chunk
  if (n > max_pairs && length(current_genes) == 0) {
    chunks[[sprintf("gene_chunk_%03d", chunk_i)]] <- g
    chunk_i <- chunk_i + 1
    cat(sprintf("[warn] gene %s alone has %d pairs, exceeds max_pairs=%d, placed in its own chunk\n",
                g, n, max_pairs))
    next
  }

  would_exceed_pairs <- (current_pairs + n) > max_pairs
  would_exceed_genes <- length(current_genes) >= genes_per_chunk

  if ((would_exceed_pairs || would_exceed_genes) && length(current_genes) > 0) {
    chunks[[sprintf("gene_chunk_%03d", chunk_i)]] <- current_genes
    chunk_i <- chunk_i + 1
    current_genes <- character(0)
    current_pairs <- 0
  }

  current_genes <- c(current_genes, g)
  current_pairs <- current_pairs + n
}

if (length(current_genes) > 0) {
  chunks[[sprintf("gene_chunk_%03d", chunk_i)]] <- current_genes
}

for (chunk_name in names(chunks)) {
  write(chunks[[chunk_name]], file = file.path(chunk_dir, chunk_name))
}

sizes <- sapply(chunks, length)
cat(sprintf("[%s] %d genes into %d chunks | genes/chunk: min=%d median=%.0f max=%d (cap=%d genes / %d pairs)\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sum(sizes), length(chunks),
            min(sizes), median(sizes), max(sizes), genes_per_chunk, max_pairs))

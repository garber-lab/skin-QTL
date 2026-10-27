# decide gene chunks basesd on SNP:gene pairs
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(magrittr)
})

dir <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models"

# read file 
f <- fread(paste0(dir,"/data/SNPs_near_TSS.bed"))

# summarize number of SNPs per gene
summary <- f[, .(nSNP = .N), by = gene_name]
setnames(summary, "gene_name", "gene")
rm(f)

set.seed(42)
summary <- summary[sample(.N)]

genes_per_chunk <- 50
max_pairs <- 100000

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
  write(chunks[[chunk_name]], file = paste0(dir, "/data/chunks/", chunk_name))
}


sizes <- sapply(chunks, length)
cat(sprintf("[%s] %d genes into %d chunks | genes/chunk: min=%d median=%.0f max=%d\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sum(sizes), length(chunks),
            min(sizes), median(sizes), max(sizes)))


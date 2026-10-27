#!/usr/bin/env Rscript
# =============================================================================
# step9.1b_make_alt_fasta.R   (analysis #3, part 1b)
# Read the REF windows fasta (from bedtools getfasta -name), substitute the
# center base with the ALT allele to make the ALT windows fasta. Sanity-checks
# that the REF window's center base matches the recorded REF allele.
# Usage: Rscript step9.1b_make_alt_fasta.R <ref.fa> <signal_info.txt> <W> <alt.fa>
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
a <- commandArgs(trailingOnly=TRUE); reffa<-a[1]; infof<-a[2]; W<-as.integer(a[3]); altfa<-a[4]
info <- fread(infof); setkey(info, snp)
L <- readLines(reffa)
hdr <- grepl("^>", L); ids <- sub("^>","", sub("::.*","", L[hdr]))   # bedtools -name -> ">SNP::chr:start-end"
seqs <- L[!hdr]
ctr <- W + 1L
out <- character(0); nmis <- 0
for (i in seq_along(ids)) {
  s <- toupper(seqs[i]); snp <- ids[i]
  ref <- info[snp, REF]; alt <- info[snp, ALT]
  if (is.na(ref) || nchar(s) < ctr) next
  if (substr(s, ctr, ctr) != ref) { nmis <- nmis + 1 }   # strand/indel mismatch -> keep window but note
  substr(s, ctr, ctr) <- alt
  out <- c(out, paste0(">", snp), s)
}
writeLines(out, altfa)
cat(sprintf("ALT fasta: %d sequences | center REF-allele mismatches: %d (check strand if high)\n", length(ids), nmis))

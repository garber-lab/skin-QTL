#!/usr/bin/env Rscript
# =============================================================================
# step9.1_build_signal_windows.R   (analysis #3, part 1)
#
# Build +/-W bp windows around every independent cis SIGNAL SNP (step8.1) for
# allelic FIMO motif-disruption scanning. Emits a BED of windows (name = SNP) and
# a signal-info table (SNP -> genes/contexts/ranks it is a signal for, REF/ALT).
# SNVs only (indels flagged and skipped). REF/ALT come from the genotype table.
#
# Inputs : ROOT/results/conditional_analysis/within_context/signals_*.txt
#          ROOT/results/genotype_of_eGene_tested_variants.txt (CHROM POS ID REF ALT)
# Outputs: ROOT/results/conditional_analysis/fimo/signal_windows.bed
#          ROOT/results/conditional_analysis/fimo/signal_info.txt
# Usage  : Rscript step9.1_build_signal_windows.R [ROOT] [W]
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
args <- commandArgs(trailingOnly=TRUE)
ROOT <- if (length(args)>=1) args[1] else "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models"
W    <- if (length(args)>=2) as.integer(args[2]) else 30L
RES  <- file.path(ROOT,"results"); CA <- file.path(RES,"conditional_analysis")
OUT  <- file.path(CA,"fimo"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

S <- rbindlist(lapply(list.files(file.path(CA,"within_context"), full.names=TRUE), function(f){
  d<-fread(f); if(!nrow(d)) return(NULL); d[, ctx:=sub("signals_","",sub(".txt","",basename(f)))] }))
# per-SNP roll-up: which gene:ctx:rank it is a signal for
snp_roll <- S[, .(members=paste(sprintf("%s:%s:sig%d", gene, ctx, rank), collapse=";"),
                  genes=paste(sort(unique(gene)),collapse=","),
                  n_ge2gene=sum(gene %in% S[, .(k=max(rank)),.(ctx,gene)][k>=2, gene])),
              by=snp]

# positions + alleles
need <- unique(S$snp); idf<-tempfile(); writeLines(need, idf)
g <- fread(cmd=sprintf("awk -F'\t' 'NR==FNR{a[$1];next} FNR==1||($3 in a)' %s %s",
           shQuote(idf), shQuote(file.path(RES,"genotype_of_eGene_tested_variants.txt"))),
           select=c("CHROM","POS","ID","REF","ALT")); unlink(idf)
setnames(g, c("CHROM","POS","ID","REF","ALT"))
info <- merge(snp_roll, g, by.x="snp", by.y="ID")
info[, is_snv := nchar(REF)==1 & nchar(ALT)==1]
snv <- info[is_snv==TRUE]
# BED (0-based): window [POS-1-W, POS+W); center base = REF at 1-based index W+1
bed <- snv[, .(chr=CHROM, start=POS-1-W, end=POS+W, name=snp, score=0, strand="+")]
bed <- bed[start>=0]
fwrite(bed, file.path(OUT,"signal_windows.bed"), sep="\t", col.names=FALSE)
fwrite(info, file.path(OUT,"signal_info.txt"), sep="\t")
cat(sprintf("signals: %d unique SNPs | SNVs windowed: %d | indels skipped: %d | center index = %d\n",
            nrow(info), nrow(bed), info[is_snv==FALSE,.N], W+1L))

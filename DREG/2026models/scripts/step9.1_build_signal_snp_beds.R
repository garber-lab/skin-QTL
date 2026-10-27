#!/usr/bin/env Rscript
# =============================================================================
# step9.1_build_signal_snp_beds.R   (analysis #3, part 1)
#
# For every independent cis SIGNAL SNP of a >=2-signal gene (step8.1), write a
# per-SNP BED (chrom start end id REF ALT; 0-based start) that the proven per-SNP
# FIMO shell driver consumes. Also writes the SNP list and a snp->gene:ctx:rank map.
#
# Inputs : ROOT/results/conditional_analysis/within_context/signals_*.txt
#          ROOT/results/genotype_of_eGene_tested_variants.txt (CHROM POS ID REF ALT)
# Outputs: ROOT/results/conditional_analysis/fimo/beds/QTL_<snp>.bed  (one per SNP)
#          ROOT/results/conditional_analysis/fimo/signal_snp_list.txt
#          ROOT/results/conditional_analysis/fimo/signal_info.txt
# Usage  : Rscript step9.1_build_signal_snp_beds.R [ROOT]
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
args <- commandArgs(trailingOnly=TRUE)
ROOT <- if (length(args)>=1) args[1] else "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models"
RES  <- file.path(ROOT,"results"); CA <- file.path(RES,"conditional_analysis")
FIMODIR <- file.path(CA,"fimo"); BEDS <- file.path(FIMODIR,"beds")
dir.create(BEDS, showWarnings=FALSE, recursive=TRUE)

S <- rbindlist(lapply(list.files(file.path(CA,"within_context"), full.names=TRUE), function(f){
  d<-fread(f); if(!nrow(d)) return(NULL); d[, ctx:=sub("signals_","",sub(".txt","",basename(f)))] }))
S[, k := max(rank), by=.(ctx,gene)]
multi <- S[k>=2]                                  # signals of >=2-signal genes only
snps  <- unique(multi$snp)

# snp -> which gene:ctx:rank it is a signal for (for step9.2 co-signal grouping)
info <- multi[, .(members=paste(sprintf("%s:%s:sig%d", gene, ctx, rank), collapse=";"),
                  genes=paste(sort(unique(gene)),collapse=",")), by=snp]

# positions + alleles
idf<-tempfile(); writeLines(snps, idf)
g <- fread(cmd=sprintf("awk -F'\t' 'NR==FNR{a[$1];next} FNR==1||($3 in a)' %s %s",
           shQuote(idf), shQuote(file.path(RES,"genotype_of_eGene_tested_variants.txt"))),
           select=c("CHROM","POS","ID","REF","ALT")); unlink(idf)
g <- g[ID %in% snps]
info <- merge(info, g, by.x="snp", by.y="ID")

# one BED per SNP: chrom  POS-1(0-based)  POS  snp  REF  ALT
for (i in seq_len(nrow(info))) {
  b <- info[i]
  fwrite(data.table(b$CHROM, b$POS-1L, b$POS, b$snp, b$REF, b$ALT),
         file.path(BEDS, sprintf("QTL_%s.bed", b$snp)), sep="\t", col.names=FALSE)
}
writeLines(info$snp, file.path(FIMODIR,"signal_snp_list.txt"))
fwrite(info, file.path(FIMODIR,"signal_info.txt"), sep="\t")
cat(sprintf(">=2-signal genes: %d | signal SNPs: %d | BEDs written: %d\n",
            uniqueN(paste(multi$ctx,multi$gene)), length(snps), nrow(info)))

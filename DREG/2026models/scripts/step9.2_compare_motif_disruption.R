#!/usr/bin/env Rscript
# =============================================================================
# step9.2_compare_motif_disruption.R   (analysis #3, part 2)
#
# Aggregate the per-SNP allelic FIMO comparisons (step9.1 -> fimo/out/
# fimo_output_<snp>.txt) and, for genes with >=2 independent signals, ask whether
# the co-signals disrupt DIFFERENT TFs (division of labour) or the SAME TF
# (reinforcement); tabulate co-disrupted TF pairs across co-signals (synergy).
#
# "Disrupted" motif at a SNP = motif_change_class in {lost_in_ALT, gained_in_ALT}
# OR |delta_log10p_ALT_minus_REF| >= DELTA (default 2 = 100x p change).
# TF symbol = prefix of the HOCOMOCO motif_id (e.g. IRF1.H12.0.A -> IRF1).
#
# Inputs : ROOT/results/conditional_analysis/fimo/out/fimo_output_<snp>.txt
#          ROOT/results/conditional_analysis/within_context/signals_*.txt
# Outputs: fimo/signal_motif_disruption.txt   (per SNP x TF, disrupted only)
#          fimo/cosignal_TF_comparison.txt      (per >=2-signal gene: verdict, TFs)
#          fimo/cosignal_TF_pairs.txt           (co-disrupted TF pairs -> synergy)
# Usage  : Rscript step9.2_compare_motif_disruption.R [ROOT] [DELTA]
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
args  <- commandArgs(trailingOnly=TRUE)
ROOT  <- if (length(args)>=1) args[1] else "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models"
DELTA <- if (length(args)>=2) as.numeric(args[2]) else 2
CA <- file.path(ROOT,"results","conditional_analysis"); FIMODIR <- file.path(CA,"fimo")

ff <- list.files(file.path(FIMODIR,"out"), pattern="^fimo_output_.*\\.txt$", full.names=TRUE)
if (!length(ff)) stop("no fimo_output_*.txt in ", file.path(FIMODIR,"out"))
D <- rbindlist(lapply(ff, function(f){ x<-fread(f); if(!nrow(x)) return(NULL)
  x[, snp := sub("^fimo_output_","",sub("\\.txt$","",basename(f)))]; x }), fill=TRUE)
D[, TF := sub("\\..*","", motif_id)]
D[, disrupted := motif_change_class %in% c("lost_in_ALT","gained_in_ALT") |
     (is.finite(delta_log10p_ALT_minus_REF) & abs(delta_log10p_ALT_minus_REF) >= DELTA)]
fwrite(D[disrupted==TRUE, .(snp, TF, motif_id, motif_change_class,
        delta_log10p_ALT_minus_REF, delta_score_ALT_minus_REF)][order(snp,-abs(delta_log10p_ALT_minus_REF))],
       file.path(FIMODIR,"signal_motif_disruption.txt"), sep="\t")

# ---- co-signal comparison per >=2-signal gene ------------------------------
S <- rbindlist(lapply(list.files(file.path(CA,"within_context"), full.names=TRUE), function(f){
  d<-fread(f); if(!nrow(d)) return(NULL); d[, ctx:=sub("signals_","",sub(".txt","",basename(f)))] }))
S[, k := max(rank), by=.(ctx,gene)]; multi <- S[k>=2]
disr <- D[disrupted==TRUE, .(tfs=list(sort(unique(TF)))), by=snp]; setkey(disr, snp)
get_tfs <- function(x){ v<-disr[x, tfs][[1]]; if(is.null(v)) character(0) else v }

cmp <- multi[, {
  tflist <- lapply(snp, get_tfs)
  n_sig_disrupt <- sum(vapply(tflist, length, 1L) > 0)
  allp <- if(length(tflist)>=2) combn(length(tflist),2,simplify=FALSE) else list()
  shared <- unique(unlist(lapply(allp, function(p) intersect(tflist[[p[1]]], tflist[[p[2]]]))))
  uni <- unique(unlist(tflist))
  .(n_signals=.N, n_sig_disrupt=n_sig_disrupt,
    signals=paste(sprintf("sig%d:%s",rank,snp),collapse=";"),
    n_TF_total=length(uni), n_TF_shared=length(shared),
    shared_TFs=paste(shared,collapse=","), all_TFs=paste(uni,collapse=","),
    verdict=fifelse(length(uni)==0,"no_disruption",
             fifelse(n_sig_disrupt<2,"single_signal_disrupts",
              fifelse(length(shared)==length(uni),"same_TF(reinforce)",
               fifelse(length(shared)==0,"different_TF(division)","mixed")))))
}, by=.(ctx,gene)]
fwrite(cmp[order(gene)], file.path(FIMODIR,"cosignal_TF_comparison.txt"), sep="\t")

# ---- co-disrupted TF pairs across co-signals (synergy candidates) ----------
pairs <- multi[, {
  tflist <- lapply(snp, get_tfs); out <- list()
  if (length(tflist)>=2) for (p in combn(length(tflist),2,simplify=FALSE))
    for (ta in tflist[[p[1]]]) for (tb in tflist[[p[2]]]) if (ta!=tb)
      out[[length(out)+1]] <- data.table(TF_A=pmin(ta,tb), TF_B=pmax(ta,tb))
  if (length(out)) rbindlist(out) else NULL
}, by=.(ctx,gene)]
if (nrow(pairs)) {
  ptab <- unique(pairs)[, .(n_genes=uniqueN(gene), genes=paste(sort(unique(gene)),collapse=",")),
                        by=.(TF_A,TF_B)][order(-n_genes)]
  fwrite(ptab, file.path(FIMODIR,"cosignal_TF_pairs.txt"), sep="\t")
}
cat("=== co-signal verdicts (>=2-signal genes) ===\n"); print(cmp[, .N, verdict][order(-N)])
cat("\n=== top co-disrupted TF pairs (synergy candidates) ===\n"); if (exists("ptab")) print(head(ptab,20))
cat("\nwrote signal_motif_disruption.txt, cosignal_TF_comparison.txt, cosignal_TF_pairs.txt\n")

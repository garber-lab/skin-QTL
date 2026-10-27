#!/usr/bin/env Rscript
# =============================================================================
# step8.2_pairwise_signal_matching.R   (genome-wide run)
#
# STEP 2. Using the within-context signal sets from step8.1, compare pairs of
# contexts (two CONDITIONS within a celltype, or two CELLTYPES within a
# condition) and place each context's PRIMARY lead within the OTHER context's
# signal hierarchy. For query-lead Q (e.g. B's primary) vs reference context R
# (e.g. A) with signals {r1,r2,...}:
#   - LD r2(Q, r_k): if max r2 >= LD_SAME -> Q is R's signal #k
#         k==1 -> same primary ;  k>1 -> Q is a SECONDARY signal of R (rank shift)
#   - else test Q in R's data:  p_marg(Q in R)  and  p_cond(Q | {r1..rK} in R)
#         p_marg not sig            -> Q is ABSENT in R  (context-specific element)
#         p_marg sig & p_cond sig   -> Q is an INDEPENDENT signal in R (distinct causal)
#         p_marg sig & p_cond n.s.  -> Q tagged by R's signals combined (complex LD)
# Done reciprocally (B-primary vs A, and A-primary vs B). Raw stats out; the
# BH-corrected final call + summary is step8.3.
#
# Inputs (cluster): VST, metadata, per-ct exprPCs, genotype PCs, genotype matrix;
#   results/conditional_analysis/within_context/signals_{ct}_{cond}.txt (step8.1)
#   results/{ct}_{cond}_eQTL_gene_fdr05_genelist.txt (union pool)
# Output: results/conditional_analysis/pairwise/matched_{ctxA}_vs_{ctxB}.txt
#
# Usage: Rscript step8.2_pairwise_signal_matching.R [ROOT] [nExprPC] [nGPC] [LD_SAME]
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
stamp <- function(...) message(sprintf("[%s] %s", format(Sys.time(),"%Y-%m-%d %H:%M:%S"), paste0(...)))
`%||%` <- function(a,b) if (is.null(a)||length(a)==0||is.na(a[1])) b else a

args    <- commandArgs(trailingOnly = TRUE)
ROOT    <- if (length(args)>=1) args[1] else "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models"
nExprPC <- if (length(args)>=2) as.integer(args[2]) else 10L
nGPC    <- if (length(args)>=3) as.integer(args[3]) else 3L
LD_SAME <- if (length(args)>=4) as.numeric(args[4]) else 0.8
QTL <- "eQTL"; celltypes <- c("FRB","KRT","MEL"); conds <- c("PBS","IFNG","IFNB","TNF")

RESULTS <- file.path(ROOT,"results")
SIGDIR  <- file.path(RESULTS,"conditional_analysis","within_context")
OUT     <- file.path(RESULTS,"conditional_analysis","pairwise"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
VST_FILE      <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
META_FILE     <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
EXPR_PC_FILE  <- function(ct) sprintf("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_opt/covariates/expression_pcs_condremoved/%s/exprPC_K%02d.tsv", ct, nExprPC)
GENO_PCS_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/ancestry_PCA_rerun2026/output/genotype_PCs_for_modeling_imputed.txt"
GENO_FILE     <- file.path(RESULTS,"genotype_of_eGene_tested_variants.txt")
ctx <- function(ct,co) paste(ct,co,sep="_")

## ---- load data + signals ---------------------------------------------------
stamp("loading data + step8.1 signals ...")
meta <- fread(META_FILE)[celltype %in% celltypes, .(sample,donor,condition,celltype)]
egenes <- list(); signals <- list()
for (ct in celltypes) for (co in conds) {
  k <- ctx(ct,co)
  gl <- file.path(RESULTS, sprintf("%s_%s_%s_gene_fdr05_genelist.txt", ct, co, QTL))
  egenes[[k]] <- if (file.exists(gl)) unique(trimws(readLines(gl))) else character(0)
  sf <- file.path(SIGDIR, sprintf("signals_%s.txt", k))
  signals[[k]] <- if (file.exists(sf)) fread(sf) else data.table()
}
union_genes <- sort(unique(unlist(egenes)))
vst <- fread(VST_FILE); setnames(vst,1,"gene"); vst <- vst[gene %in% union_genes]
gcol <- setdiff(names(vst),"gene"); VM <- as.matrix(vst[,..gcol]); rownames(VM) <- vst$gene
exprpc <- list(); epc_cols <- list()
for (ct in celltypes){ e <- fread(EXPR_PC_FILE(ct)); setnames(e,1,"sample"); setkey(e,sample); exprpc[[ct]] <- e; epc_cols[[ct]] <- setdiff(names(e),"sample") }
gp <- fread(GENO_PCS_FILE); gpc_cols <- paste0("PC",seq_len(nGPC)); gp <- gp[,c("donor",gpc_cols),with=FALSE]; setkey(gp,donor)
stamp("loading genotype matrix ...")
G <- fread(GENO_FILE); dcol <- setdiff(names(G),c("CHROM","POS","ID","REF","ALT"))
GM <- as.matrix(G[,..dcol]); rownames(GM) <- G$ID; mode(GM) <- "numeric"; rm(G); gc()

ctx_samples <- list()
for (ct in celltypes) for (co in conds) {
  s <- meta[celltype==ct & condition==co, .(sample,donor)]
  ctx_samples[[ctx(ct,co)]] <- s[sample %in% colnames(VM) & donor %in% colnames(GM) & donor %in% gp$donor & sample %in% exprpc[[ct]]$sample]
}
ld2 <- function(a,b){ if(is.na(a)||is.na(b)||!(a %in% rownames(GM))||!(b %in% rownames(GM))) return(NA_real_)
  ga<-GM[a,]; gb<-GM[b,]; ok<-is.finite(ga)&is.finite(gb); if(sum(ok)<3||sd(ga[ok])==0||sd(gb[ok])==0) return(NA_real_); cor(ga[ok],gb[ok])^2 }

# marginal & conditional (on cond_snps set) p of test_snp for gene in context ct,co
assoc <- function(gene, test_snp, cond_snps, ct, co) {
  s <- ctx_samples[[ctx(ct,co)]]
  if (nrow(s)<10 || !(gene %in% rownames(VM)) || is.na(test_snp) || !(test_snp %in% rownames(GM))) return(list(p_marg=NA_real_,p_cond=NA_real_))
  d <- data.table(sample=s$sample, donor=s$donor, vst=as.numeric(VM[gene,s$sample]), g_test=GM[test_snp,s$donor])
  d <- merge(merge(d, exprpc[[ct]], by="sample"), gp, by="donor")
  covs <- c(epc_cols[[ct]], gpc_cols)
  d <- d[complete.cases(d[,c("vst","g_test",covs),with=FALSE])]
  if (length(unique(d$g_test))<2 || nrow(d)<length(covs)+6) return(list(p_marg=NA_real_,p_cond=NA_real_))
  f0 <- lm(reformulate(c("g_test",covs),"vst"), d); pm <- coef(summary(f0))["g_test","Pr(>|t|)"]
  cs <- setdiff(intersect(cond_snps, rownames(GM)), test_snp)
  pc <- NA_real_
  if (length(cs)) {
    cg <- t(GM[cs, d$donor, drop=FALSE]); colnames(cg) <- paste0("C",seq_along(cs))
    dd <- cbind(d, as.data.table(cg)); cgcols <- colnames(cg)
    dd <- dd[complete.cases(dd[,c("vst","g_test",covs,cgcols),with=FALSE])]
    keepc <- cgcols[sapply(cgcols, function(cc) length(unique(dd[[cc]]))>=2)]
    if (nrow(dd)>=length(covs)+length(keepc)+6 && length(unique(dd$g_test))>=2) {
      f1 <- lm(reformulate(c("g_test",keepc,covs),"vst"), dd); pc <- coef(summary(f1))["g_test","Pr(>|t|)"]
    }
  }
  list(p_marg=pm, p_cond=pc)
}

# classify query lead vs reference context's signal set
classify <- function(gene, query, ref_ct, ref_co, ref_sig) {
  if (is.na(query)) return(data.table(match_rank=NA_integer_, best_r2=NA_real_, p_marg=NA_real_, p_cond=NA_real_, class_pre="no_query_signal"))
  r2 <- if (nrow(ref_sig)) vapply(ref_sig$snp, function(s) ld2(query,s), numeric(1)) else numeric(0)
  best_r2 <- if(length(r2)) max(r2, na.rm=TRUE) else NA_real_
  best_rank <- if(length(r2) && is.finite(best_r2)) ref_sig$rank[which.max(r2)] else NA_integer_
  if (is.finite(best_r2) && best_r2 >= LD_SAME) {
    cls <- if (best_rank==1) "same_primary" else "matches_ref_secondary"
    return(data.table(match_rank=best_rank, best_r2=round(best_r2,3), p_marg=NA_real_, p_cond=NA_real_, class_pre=cls))
  }
  a <- assoc(gene, query, ref_sig$snp, ref_ct, ref_co)
  data.table(match_rank=best_rank, best_r2=round(best_r2,3), p_marg=a$p_marg, p_cond=a$p_cond, class_pre="test_needed")
}

run_pair <- function(ctA,coA,ctB,coB){
  kA<-ctx(ctA,coA); kB<-ctx(ctB,coB)
  mode <- if (ctA==ctB) "across_conditions" else "across_celltypes"; held <- if (ctA==ctB) ctA else coA
  pool <- union(egenes[[kA]], egenes[[kB]])
  sigA <- signals[[kA]]; sigB <- signals[[kB]]
  rows <- vector("list", length(pool))
  for (i in seq_along(pool)) {
    g <- pool[i]
    aSig <- if (nrow(sigA)) sigA[gene==g][order(rank)] else data.table(); bSig <- if (nrow(sigB)) sigB[gene==g][order(rank)] else data.table()
    a1 <- if (nrow(aSig)) aSig$snp[1] else NA_character_
    b1 <- if (nrow(bSig)) bSig$snp[1] else NA_character_
    if (is.na(a1) && is.na(b1)) next
    cB <- classify(g, b1, ctA, coA, aSig)   # B's primary vs A
    cA <- classify(g, a1, ctB, coB, bSig)   # A's primary vs B
    rows[[i]] <- data.table(mode, held, ctxA=kA, ctxB=kB, gene=g,
      a1=a1, nA_signals=nrow(aSig), b1=b1, nB_signals=nrow(bSig), same_lead=(!is.na(a1)&&!is.na(b1)&&a1==b1),
      b1vsA_class=cB$class_pre, b1vsA_rank=cB$match_rank, b1vsA_r2=cB$best_r2, b1vsA_pmarg=cB$p_marg, b1vsA_pcond=cB$p_cond,
      a1vsB_class=cA$class_pre, a1vsB_rank=cA$match_rank, a1vsB_r2=cA$best_r2, a1vsB_pmarg=cA$p_marg, a1vsB_pcond=cA$p_cond)
  }
  out <- rbindlist(rows)
  fwrite(out, file.path(OUT, sprintf("matched_%s_vs_%s.txt", kA, kB)), sep="\t")
  stamp("  ",mode,"  ",kA," vs ",kB,": ",nrow(out)," genes")
}

stamp("running pairwise matching ...")
for (ct in celltypes){ cp<-combn(conds,2); for(j in seq_len(ncol(cp))) run_pair(ct,cp[1,j],ct,cp[2,j]) }
for (co in conds){ cp<-combn(celltypes,2); for(j in seq_len(ncol(cp))) run_pair(cp[1,j],co,cp[2,j],co) }
stamp("DONE step8.2.")

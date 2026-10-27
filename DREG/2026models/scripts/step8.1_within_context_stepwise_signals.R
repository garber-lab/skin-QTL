#!/usr/bin/env Rscript
# =============================================================================
# step8.1_within_context_stepwise_signals.R   (genome-wide run)
#
# STEP 1 of the multi-signal conditional pipeline. Within each context
# (celltype x condition), for each eGene, discover the INDEPENDENT cis signals by
# forward stepwise conditional analysis (GCTA-COJO style, exact with individual
# genotypes): signal-1 = lead; condition on it and re-scan all cis SNPs for
# signal-2; condition on 1+2 for signal-3; ... until no SNP passes the per-gene
# threshold (or maxK reached). Fast via matrix residualization (qr.resid).
#
# Model covariates match step1: intercept + exprPC1..K(condition-removed, per
# celltype) + genotypePC1..nGPC. One sample per donor per context -> lm/OLS.
# Acceptance threshold per gene: Bonferroni alpha / n_cis (default alpha=0.05).
#
# Inputs (cluster): VST, metadata, per-ct exprPCs, genotype PCs,
#   results/{ct}_{cond}_eQTL_gene_fdr05_genelist.txt   (eGenes to scan)
#   results/QTL_output_of_all_tested_variants_of_union_of_eGenes_{ct}_{cond}_eQTL.txt  (cis SNP list per gene)
#   results/genotype_of_eGene_tested_variants.txt      (donor genotypes)
# Output: results/conditional_analysis/within_context/signals_{ct}_{cond}.txt
#           (gene, rank, snp, dist, cond_p, beta, n_cis, n_donors)
#
# Usage: Rscript step8.1_within_context_stepwise_signals.R [ROOT] [nExprPC] [nGPC] [maxK] [alpha]
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
stamp <- function(...) message(sprintf("[%s] %s", format(Sys.time(),"%Y-%m-%d %H:%M:%S"), paste0(...)))

args    <- commandArgs(trailingOnly = TRUE)
ROOT    <- if (length(args)>=1) args[1] else "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models"
nExprPC <- if (length(args)>=2) as.integer(args[2]) else 10L
nGPC    <- if (length(args)>=3) as.integer(args[3]) else 3L
maxK    <- if (length(args)>=4) as.integer(args[4]) else 5L
alpha   <- if (length(args)>=5) as.numeric(args[5]) else 0.05
QTL <- "eQTL"; celltypes <- c("FRB","KRT","MEL"); conds <- c("PBS","IFNG","IFNB","TNF")

RESULTS <- file.path(ROOT,"results")
OUT <- file.path(RESULTS,"conditional_analysis","within_context"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
VST_FILE      <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
META_FILE     <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
EXPR_PC_FILE  <- function(ct) sprintf("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_opt/covariates/expression_pcs_condremoved/%s/exprPC_K%02d.tsv", ct, nExprPC)
GENO_PCS_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/ancestry_PCA_rerun2026/output/genotype_PCs_for_modeling_imputed.txt"
GENO_FILE     <- file.path(RESULTS,"genotype_of_eGene_tested_variants.txt")

## ---- load shared data ------------------------------------------------------
stamp("loading metadata / VST / covariates / genotypes ...")
meta <- fread(META_FILE)[celltype %in% celltypes, .(sample,donor,condition,celltype)]
egenes <- list()
for (ct in celltypes) for (co in conds) {
  gl <- file.path(RESULTS, sprintf("%s_%s_%s_gene_fdr05_genelist.txt", ct, co, QTL))
  egenes[[paste(ct,co,sep="_")]] <- if (file.exists(gl)) unique(trimws(readLines(gl))) else character(0)
}
union_genes <- sort(unique(unlist(egenes)))
vst <- fread(VST_FILE); setnames(vst,1,"gene"); vst <- vst[gene %in% union_genes]
gcol <- setdiff(names(vst),"gene"); VM <- as.matrix(vst[,..gcol]); rownames(VM) <- vst$gene
exprpc <- list(); epc_cols <- list()
for (ct in celltypes){ e <- fread(EXPR_PC_FILE(ct)); setnames(e,1,"sample"); setkey(e,sample); exprpc[[ct]] <- e; epc_cols[[ct]] <- setdiff(names(e),"sample") }
gp <- fread(GENO_PCS_FILE); gpc_cols <- paste0("PC",seq_len(nGPC)); gp <- gp[,c("donor",gpc_cols),with=FALSE]; setkey(gp,donor)
stamp("loading genotype matrix (~1.2GB) ...")
G <- fread(GENO_FILE); dcol <- setdiff(names(G),c("CHROM","POS","ID","REF","ALT"))
GM <- as.matrix(G[,..dcol]); rownames(GM) <- G$ID; mode(GM) <- "numeric"; rm(G); gc()
stamp("genotype matrix: ", nrow(GM)," x ",ncol(GM))

## ---- stepwise conditional signal finder ------------------------------------
# y (n), X (n x p, incl intercept), Gm (n x m cis dosages, mean-imputed), returns signals
stepwise <- function(y, X, Gm, snp_ids, dist_v) {
  m <- ncol(Gm); if (m < 1) return(NULL)
  alpha_g <- alpha / m                                   # per-gene Bonferroni
  sel <- integer(0); sigs <- list()
  for (k in seq_len(maxK)) {
    Xk <- if (length(sel)) cbind(X, Gm[, sel, drop=FALSE]) else X
    q  <- qr(Xk)
    ry <- qr.resid(q, y); rG <- qr.resid(q, Gm)
    df <- length(y) - ncol(Xk) - 1; if (df < 5) break
    sg <- colSums(rG^2); ss_y <- sum(ry^2)
    b  <- as.numeric(crossprod(rG, ry)) / sg
    rss <- ss_y - b^2*sg
    se  <- sqrt((rss/df)/sg)
    p   <- 2*pt(-abs(b/se), df)
    p[sel] <- Inf; p[!is.finite(p) | sg < 1e-8] <- Inf
    j <- which.min(p)
    if (!length(j) || p[j] >= alpha_g) break
    sigs[[k]] <- data.table(rank=k, snp=snp_ids[j], dist=dist_v[j], cond_p=p[j], beta=b[j])
    sel <- c(sel, j)
  }
  if (length(sigs)) rbindlist(sigs) else NULL
}

impute_cols <- function(M) { cm <- colMeans(M, na.rm=TRUE); i <- which(is.na(M), arr.ind=TRUE)
  if (nrow(i)) M[i] <- cm[i[,2]]; M }

## ---- process context by context --------------------------------------------
for (ct in celltypes) for (co in conds) {
  key <- paste(ct,co,sep="_"); eg <- egenes[[key]]; if (!length(eg)) next
  uf <- file.path(RESULTS, sprintf("QTL_output_of_all_tested_variants_of_union_of_eGenes_%s_%s_%s.txt", ct, co, QTL))
  if (!file.exists(uf)) { stamp("  [",key,"] missing union file -- skip"); next }
  cis <- fread(uf, select=c("snp","gene","dist"))[gene %in% eg & abs(dist)<=5e5]
  cis_by <- split(cis, cis$gene)
  s <- meta[celltype==ct & condition==co, .(sample,donor)]
  s <- s[sample %in% colnames(VM) & donor %in% colnames(GM) & donor %in% gp$donor & sample %in% exprpc[[ct]]$sample]
  epc_df <- as.data.frame(exprpc[[ct]]); gp_df <- as.data.frame(gp)
  Emat <- as.matrix(epc_df[match(s$sample, epc_df$sample), epc_cols[[ct]], drop=FALSE])
  Pmat <- as.matrix(gp_df[match(s$donor, gp_df$donor), gpc_cols, drop=FALSE])
  Xcov <- cbind(1, Emat, Pmat)
  res <- list()
  for (g in eg) {
    cc <- cis_by[[g]]; if (is.null(cc)) next
    snps <- intersect(cc$snp, rownames(GM)); if (length(snps) < 1) next
    cc <- cc[match(snps, snp)]
    y <- as.numeric(VM[g, s$sample])
    Gm <- t(GM[snps, s$donor, drop=FALSE]); Gm <- impute_cols(Gm)   # donors x snps
    keep <- is.finite(y) & complete.cases(Xcov)
    if (sum(keep) < ncol(Xcov)+6) next
    sg <- stepwise(y[keep], Xcov[keep,,drop=FALSE], Gm[keep,,drop=FALSE], cc$snp, cc$dist)
    if (!is.null(sg)) { sg[, `:=`(gene=g, n_cis=length(snps), n_donors=sum(keep))]; res[[g]] <- sg }
  }
  out <- rbindlist(res, use.names=TRUE)
  if (nrow(out)) setcolorder(out, c("gene","rank","snp","dist","cond_p","beta","n_cis","n_donors"))
  fwrite(out, file.path(OUT, sprintf("signals_%s.txt", key)), sep="\t")
  stamp("  [",key,"] eGenes ",length(eg),
        " -> genes w/ signals ", if(nrow(out)) uniqueN(out$gene) else 0,
        ", total signals ", nrow(out),
        " (>=2 signals: ", if(nrow(out)) out[, .N, gene][N>=2, .N] else 0, " genes)")
}
stamp("DONE step8.1.")

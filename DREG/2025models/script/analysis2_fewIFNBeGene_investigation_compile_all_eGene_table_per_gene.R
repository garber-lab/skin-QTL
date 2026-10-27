#!/usr/bin/env Rscript
# =============================================================================
# Per-gene builder for the "why so few IFNB eGenes" analysis.
#
# For one gene, emits THREE files into working_dir:
#   leads.txt          - lead tie-BLOCKS: every variant sharing the best p in each
#                        celltype x cytokine anchor (keep ALL ties, not one rep).
#                        Carries tie_block_id, n_lead_tied, cbp, maf, dist.
#                        This is the object the enrichment test consumes.
#   lead_matrix.txt(T1)- each distinct lead variant re-evaluated across ALL 12
#                        contexts (incl PBS fill) + reQTL annotation. Mechanism /
#                        cross-context comparison (e.g. IFNG lead seen in IFNB).
#   window_summary.txt - gene x context window summary (T2). The count/attrition
#       (T2)            spine for "is this an eGene here, and why".
#
# Usage:  Rscript build_lead_tables.R <GENE> <working_dir>
# =============================================================================
Sys.setenv(TZ = "America/New_York")
suppressPackageStartupMessages(library(data.table))
ts <- function() format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")

args        <- commandArgs(trailingOnly = TRUE)
this_gene   <- args[1]
working_dir <- args[2]
# debug: 
# this_gene <- "GBP3"; working_dir <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats/G/GBP3"

infile   <- file.path(working_dir, "big_table.txt")
out_blk  <- file.path(working_dir, "leads.txt")
out_lead <- file.path(working_dir, "lead_matrix.txt")
out_win  <- file.path(working_dir, "window_summary.txt")

# already done? skip
if (all(file.exists(out_blk, out_lead, out_win))) {          # idempotent reruns
  cat(ts(), this_gene, "outputs exist, skipping\n"); quit(save = "no")
}

# no input? skip
if (!file.exists(infile)) {
  cat(ts(), this_gene, "big_table.txt doesn't exist, skipping\n")
  quit(save = "no")
}

# ---- knobs -----------------------------------------------------------------
CYTOKINES            <- c("IFNB", "IFNG", "TNF")  # anchors; PBS is fill-only
ANCHOR_ON_EGENE_ONLY <- TRUE                      # FALSE = anchor on every cytokine context
PTOL                 <- 1e-6                      # rel. tol for "same best p" (perfect-LD proxies)
CBP                  <- "abs_logfc_x_active_allele_quantile.mean.pval"  # primary CBP metric (lower=stronger)

safe_min      <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else min(x) }
maf_from_geno <- function(g0, g1, g2) { n <- g0 + g1 + g2; af <- (g1 + 2 * g2) / (2 * n); pmin(af, 1 - af) }

big  <- fread(infile)
eqtl <- big[QTLtype == "eQTL"]                    # eGene count = eQTL only; reQTL is annotation

# ---- 1. lead tie-block per celltype x cytokine (keep ALL tied-best) --------
blocks <- list()
for (ct in c("FRB", "KRT", "MEL")) for (cond in CYTOKINES) {
  sub <- eqtl[celltype == ct & condition == cond]
  if (!nrow(sub)) next
  if (ANCHOR_ON_EGENE_ONLY && !any(sub$is_egene %in% TRUE)) next
  pmin <- min(sub$p, na.rm = TRUE)
  tie  <- sub[p <= pmin * (1 + PTOL)]  # all perfect-LD proxies, not one representative. p-values are stores as floating-point decimals, and perfect LD variants that should be identical can differ int he 15th digit due to rounding. p == pmin would silently drop some of them.
  blocks[[paste(ct, cond)]] <- data.table(
    gene            = this_gene,
    celltype        = ct,
    condition       = cond,
    tie_block_id    = paste(this_gene, ct, cond, sep = "_"),
    anchor_context  = paste(ct, cond, "eQTL", sep = "_"),
    snp             = tie$snp,
    p               = tie$p,
    beta            = tie$beta,
    dist            = tie$dist,
    maf             = maf_from_geno(tie$n_geno0, tie$n_geno1, tie$n_geno2),
    cbp_pval        = tie[[CBP]],     # copying the column for CBP as cbp_pval, NA for IFNB/TNF (CBP only PBS/IFNG)
    ATACpeak_isopen = tie$ATACpeak_isopen,
    n_lead_tied     = nrow(tie))
}
blocks <- rbindlist(blocks)
if (!nrow(blocks)) { cat(ts(), this_gene, "no eGene cytokine context\n"); quit(save = "no") }
fwrite(blocks, out_blk, sep = "\t")

# distinct lead variants + which contexts they lead
# group blocks by snp; for each distinct variant, build two summary columns (what contexts is this variant lead in, and how many variants are tied with this one)
lead_info <- blocks[, .(lead_for_contexts = paste(sort(unique(anchor_context)), collapse = ";"),
                        n_lead_tied        = max(n_lead_tied)), by = snp]

# ---- 2. T1: each lead variant across ALL contexts (incl PBS) ---------------
T1 <- merge(eqtl[snp %in% lead_info$snp], lead_info, by = "snp", all.x = TRUE) # filter by lead variants only, then left join with lead_info
T1[, this_context    := paste(celltype, condition, "eQTL", sep = "_")]
T1[, is_context_lead := mapply(function(c, s) c %in% strsplit(s, ";")[[1]],
                               this_context, lead_for_contexts)]
T1[, gene := this_gene]

# reQTL annotation -> constitutive (eQTL in PBS, no reQTL) vs response (reQTL) axis
reqtl <- unique(big[QTLtype == "reQTL",
                    .(snp, celltype, condition, reQTL_p = p, reQTL_beta = beta, is_reQTL_sig = is_sig_snp)])
T1 <- merge(T1, reqtl, by = c("snp", "celltype", "condition"), all.x = TRUE)
fwrite(T1, out_lead, sep = "\t")

# ---- 3. T2: gene x context window summary (all 12 contexts) ----------------
win <- eqtl[, .(
  is_egene               = any(is_egene %in% TRUE),
  best_p_in_window       = safe_min(p),
  n_snps_passing_nominal = sum(is_sig_snp %in% TRUE),
  best_chrombpnet_pval   = safe_min(get(CBP)),      # lower = stronger; NA in IFNB/TNF
  any_open_ATAC          = any(ATACpeak_isopen %in% TRUE),
  any_dynamic_ATAC       = any(ATACpeak_padj < 0.05, na.rm = TRUE),
  meanCPM                = meanCPM[1],
  gene_DESeq2_log2FC     = gene_DESeq2_log2FC[1],
  gene_DESeq2_padj       = gene_DESeq2_padj[1],
  gene_category          = gene_category[1]
), by = .(celltype, condition)]
win[, gene := this_gene]
fwrite(win, out_win, sep = "\t")

cat(ts(), this_gene, "done:", uniqueN(lead_info$snp), "lead variants,",
    uniqueN(blocks$tie_block_id), "cytokine contexts in which", this_gene, "is an eGene, out of", nrow(win), "contexts\n")

message(paste0("output generated in: \n", 
        out_blk, "\n",
        out_lead, "\n",
        out_win))

#!/usr/bin/env Rscript
# =============================================================================
# Are IFNG lead variants enriched for ChromBPNet chromatin disruption vs
# MAF- and TSS-distance-matched non-lead tested variants?
#
# Design (see discussion):
#   * Unit of test = tie BLOCK (gene x celltype IFNG locus), NOT each variant,
#     so perfect-LD lead members are not pseudoreplicated.
#   * Block score = MEAN cbp over members  -> unbiased (no cherry-picking the max).
#     MAX cbp reported alongside, clearly labelled OPTIMISTIC (a single causal
#     high-CBP variant diluted by its LD neighbours otherwise).
#   * IFNG only: CBP exists for PBS/IFNG; IFNB/TNF have no score to test.
#   * Background = non-lead tested IFNG variants, matched on MAF decile x
#     log-TSS-distance tertile (CLAUDE.md backlog #5 negative control).
#
# Inputs:
#   args[1] master_leads.txt  = rbind of all per-gene leads.txt
#   args[2] gene_dirs.txt     = one big_table.txt-containing dir per line (background)
#
# Usage: Rscript test_lead_cbp_enrichment.R master_leads.txt gene_dirs.txt
# NOTE: first-pass matching (bin-mean background). Refine to 1:k nearest-neighbour
#       matching if the bins are sparse.
# =============================================================================
Sys.setenv(TZ = "America/New_York")
suppressPackageStartupMessages({ library(data.table); library(pROC) })
ts  <- function() format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")

args         <- commandArgs(trailingOnly = TRUE)
master_leads <- args[1]
gene_dirs    <- readLines(args[2])

# debug 
#master_leads="/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats/G/GBP3/master_leads.txt"; gene_dirs="/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats/gene_dirs.txt"

CBP <- "abs_logfc_x_active_allele_quantile.mean.pval"
maf_from_geno <- function(g0, g1, g2) { n <- g0 + g1 + g2; af <- (g1 + 2 * g2) / (2 * n); pmin(af, 1 - af) }
nlp <- function(p) -log10(pmax(p, .Machine$double.xmin))     # higher = stronger disruption

# ---- lead identities (IFNG only) -------------------------------------------
leads    <- fread(master_leads)[condition == "IFNG"]
lead_key <- unique(leads[, .(gene, celltype, snp)]); lead_key[, is_lead := TRUE]

# ---- stream per-gene big_tables: IFNG eQTL variants with CBP (lead+background)
pool <- vector("list", length(gene_dirs))
for (i in seq_along(gene_dirs)) {
  f <- file.path(gene_dirs[i], "big_table.txt"); if (!file.exists(f)) next
  d <- fread(f)[QTLtype == "eQTL" & condition == "IFNG" & is.finite(get(CBP)),
                .(gene, celltype, snp, dist,
                  maf = maf_from_geno(n_geno0, n_geno1, n_geno2),
                  cbp = get(CBP), inpeak = ATACpeak_isopen)]
  if (!nrow(d)) next
  d <- merge(d, lead_key, by = c("gene", "celltype", "snp"), all.x = TRUE)
  d[is.na(is_lead), is_lead := FALSE]
  pool[[i]] <- d
}
pool <- rbindlist(pool)
pool[, cbp_nlp := nlp(cbp)]
cat(ts(), "pooled", nrow(pool), "IFNG variants;", sum(pool$is_lead), "lead members\n")

# ---- collapse leads to ONE observation per tie block (gene x celltype) ------
lead_blocks <- pool[is_lead == TRUE, .(
  cbp_score = mean(cbp_nlp),          # unbiased block score
  cbp_max   = max(cbp_nlp),           # optimistic bound
  maf       = median(maf),
  dist      = median(abs(dist)),
  n_members = .N
), by = .(gene, celltype)]

# ---- matched background: MAF decile x log-dist tertile bins -----------------
maf_br  <- quantile(pool$maf, 0:10 / 10, na.rm = TRUE)
bin_maf  <- function(x) cut(x, unique(maf_br), include.lowest = TRUE, labels = FALSE)
bin_dist <- function(x) cut(log10(abs(x) + 1), 3, labels = FALSE)

bg <- pool[is_lead == FALSE]
bg[,          `:=`(mb = bin_maf(maf),  db = bin_dist(dist))]
lead_blocks[, `:=`(mb = bin_maf(maf),  db = bin_dist(dist))]
bg_bin <- bg[, .(bg_cbp_score = mean(cbp_nlp)), by = .(mb, db)]   # background equiv per bin

paired <- merge(lead_blocks, bg_bin, by = c("mb", "db"), all.x = TRUE)
paired <- paired[is.finite(bg_cbp_score)]

# ---- tests ------------------------------------------------------------------
w_mean <- wilcox.test(paired$cbp_score, paired$bg_cbp_score, paired = TRUE)
w_max  <- wilcox.test(paired$cbp_max,   paired$bg_cbp_score, paired = TRUE)

long <- rbind(data.table(grp = "lead",       v = paired$cbp_score),
              data.table(grp = "background", v = paired$bg_cbp_score))
auc  <- as.numeric(pROC::auc(response = long$grp, predictor = long$v,
                             levels = c("background", "lead"), direction = "<"))

cat(sprintf("%s blocks tested: %d\n", ts(), nrow(paired)))
cat(sprintf("  UNBIASED  (mean-cbp/block): Wilcoxon paired p = %.2e | AUROC = %.3f | median lead=%.2f bg=%.2f\n",
            w_mean$p.value, auc, median(paired$cbp_score), median(paired$bg_cbp_score)))
cat(sprintf("  OPTIMISTIC(max-cbp/block) : Wilcoxon paired p = %.2e | median lead=%.2f\n",
            w_max$p.value, median(paired$cbp_max)))

fwrite(paired, sub("\\.txt$", "_cbp_enrichment_blocks.txt", master_leads), sep = "\t")
cat(ts(), "wrote", sub("\\.txt$", "_cbp_enrichment_blocks.txt", master_leads), "\n")
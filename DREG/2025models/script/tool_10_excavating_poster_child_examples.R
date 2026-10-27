#!/usr/bin/env Rscript
# =====================================================================
# Poster-Child Story Excavation
# Genetic architecture of cytokine response in skin (FRB/KRT/MEL x PBS/IFNG)
#
# Inputs (TAB-delimited; handoff said space-delimited but the label/motif
#         columns contain spaces, so the true separator is TAB):
#   mechanistic_variant_table_detail_masterfile.txt  (snp x gene x ct x cond [x QTLtype])
#   fimo_output_masterfile.txt                       (... x TF)
#
# Outputs (all TAB-delimited):
#   top_candidates_per_context.txt   (3a)
#   variant_sharing_summary.txt      (3b)
#   tf_disruption_vocabulary.txt     (3c)
#   tss_distance_analysis.txt        (3d)
#   story_candidates_flagged.txt     (interpretation layer)
#
# DATA-REALITY NOTES (verified against the actual files, not assumed):
#   * Delimiter is TAB.
#   * is_egene, is_sig_snp (mech) and interesting_final, is_expressed (fimo)
#     are TRUE for every row, and priority_tier is only Tier1/Tier2 -> the
#     handoff's filters on these are no-ops here but are KEPT for robustness.
#   * ATACpeak_isopen has blank/NA rows -> use %in% TRUE (NA-safe), never ==.
#   * Only conditions PBS + IFNG exist (all 3 cell types) => max 6 contexts;
#     "stimulus-specific" therefore means IFNG-specific.
#   * 146 snp x gene x ct x cond keys appear twice (eQTL + reQTL): identical
#     chromatin metric, different beta. We DEDUPE to one row per variant x
#     context, PREFERRING the eQTL (main-effect) beta, so the FIMO join does
#     not double-count.
#   * strongest_change_class is only stronger/weaker_in_ALT and the gain/loss
#     COUNTS (n_gained_in_ALT, n_lost_in_ALT) are 0 for every row -> the
#     literal "Pioneer TF disruption" rule (gained/lost) matches nothing.
#     Reinterpreted as: pioneer TF_family AND top-quartile disruption strength.
#
# NAMING CAUTION (two near-identical names, different meanings):
#   * cbp_best_pval  -> a column already in the mechanistic file (secondary
#                       chromBPNet metric). Reported as-is in 3a.
#   * best_cbp_pval  -> computed here in 3b = min of the PRIMARY metric
#                       abs_logfc_x_active_allele_quantile.mean.pval across a
#                       variant's tested contexts. NOT the same column.
# =====================================================================

suppressMessages(library(data.table))

log <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ..., "\n")
safe_min <- function(x) { x <- x[!is.na(x)]; if (length(x)) min(x) else NA_real_ }

wd        <- "/Users/shuoshan/Downloads/poster_child_excavation"
mech_file <- file.path(wd, "mechanistic_variant_table_detail_masterfile.txt")
fimo_file <- file.path(wd, "fimo_output_masterfile.txt")

PRIMARY <- "abs_logfc_x_active_allele_quantile.mean.pval"   # primary chromBPNet metric
PIONEER_FAMS <- c("bHLH", "bZIP", "Forkhead", "Homeodomain", "CBF/NF-Y")

# ---------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------
log("Reading mechanistic masterfile ...")
mech <- fread(mech_file, sep = "\t", header = TRUE, na.strings = c("NA", ""))
log("  mechanistic:", nrow(mech), "rows x", ncol(mech), "cols")
stopifnot(ncol(mech) == 44L)

log("Reading FIMO masterfile (large) ...")
fimo <- fread(fimo_file, sep = "\t", header = TRUE, na.strings = c("NA", ""))
log("  fimo:", nrow(fimo), "rows x", ncol(fimo), "cols")
stopifnot(ncol(fimo) == 44L)

# Defensive boolean coercion (so %in% TRUE behaves whether read as logical or character)
for (b in c("is_egene", "is_sig_snp", "ATACpeak_isopen", "cbp_eqtl_concordant"))
  if (!is.logical(mech[[b]])) mech[, (b) := as.logical(get(b))]
for (b in c("interesting_final", "is_expressed", "is_DE", "is_induced"))
  if (!is.logical(fimo[[b]])) fimo[, (b) := as.logical(get(b))]

# ---------------------------------------------------------------------
# Dedupe eQTL/reQTL -> one row per variant x context, preferring eQTL beta
# ---------------------------------------------------------------------
mech[, QTLtype := factor(QTLtype, levels = c("eQTL", "reQTL"))]  # eQTL sorts first
setorderv(mech, c("snp", "gene", "celltype", "condition", "QTLtype"))
n_before <- nrow(mech)
mech_dedup <- unique(mech, by = c("snp", "gene", "celltype", "condition"))
log("  deduped eQTL/reQTL pairs:", n_before - nrow(mech_dedup), "rows dropped;",
    nrow(mech_dedup), "unique variant x context rows remain")

# Convenience context label
mech_dedup[, ctx := paste(celltype, condition, sep = "_")]

# ---------------------------------------------------------------------
# Step 1 -- candidate filtering + within-context ranking
# ---------------------------------------------------------------------
cand <- mech_dedup[is_egene %in% TRUE &
                   is_sig_snp %in% TRUE &
                   ATACpeak_isopen %in% TRUE &
                   cbp_eqtl_concordant %in% TRUE]
setorderv(cand, c("celltype", "condition", PRIMARY))
cand[, context_rank := seq_len(.N), by = .(celltype, condition)]
log("Candidates (egene & sig & ATAC-open & concordant):", nrow(cand),
    "rows;", uniqueN(cand[, .(gene, snp)]), "unique gene x snp;",
    "by context:")
print(cand[, .N, by = .(celltype, condition)][order(celltype, condition)])

# ---------------------------------------------------------------------
# Step 2 -- context-matched FIMO join
# ---------------------------------------------------------------------
fimo_f <- fimo[interesting_final %in% TRUE &
               is_expressed %in% TRUE &
               priority_tier %in% c("Tier1_high_confidence", "Tier2_moderate")]
fimo_sub <- fimo_f[, .(gene, snp = variant, celltype, condition,
                       TF, TF_family, strongest_change_class,
                       best_abs_delta_score, best_priority_abs_delta_log10p,
                       best_qvalue, motif_display_REF, motif_display_ALT,
                       priority_tier, effect_dir)]

joined <- merge(cand, fimo_sub,
                by = c("gene", "snp", "celltype", "condition"),
                all = FALSE)                                   # inner join (nomatch=0)
joined[, ctx := paste(celltype, condition, sep = "_")]
log("Joined candidate x FIMO-TF rows:", nrow(joined),
    "covering", uniqueN(joined[, .(gene, snp)]), "unique gene x snp")
if (nrow(joined) == 0L) stop("No candidate variants have FIMO hits -- nothing to excavate.")

# =====================================================================
# 3a -- Per-context top candidate table (best FIMO hit per variant, top 20/context)
# =====================================================================
# Deterministic best-hit-per-variant-per-context: strongest disruption first.
setorderv(joined,
          c("gene", "snp", "celltype", "condition",
            "best_priority_abs_delta_log10p", "best_abs_delta_score", "TF"),
          order = c(1, 1, 1, 1, -1, -1, 1))
best_fimo <- joined[, .SD[1L], by = .(gene, snp, celltype, condition)]

setorderv(best_fimo, c("celltype", "condition", PRIMARY, "best_priority_abs_delta_log10p"),
          order = c(1, 1, 1, -1))
top20 <- best_fimo[, head(.SD, 20L), by = .(celltype, condition)]

a_cols <- c("gene", "snp", "celltype", "condition", "beta", "dist",
            PRIMARY, "cbp_best_pval", "logfc_alt",
            "TF", "TF_family", "strongest_change_class",
            "best_abs_delta_score", "best_priority_abs_delta_log10p",
            "motif_display_REF", "motif_display_ALT", "priority_tier")
top20_out <- top20[, ..a_cols]
fwrite(top20_out, file.path(wd, "top_candidates_per_context.txt"), sep = "\t")
log("3a top_candidates_per_context.txt:", nrow(top20_out), "rows")

# =====================================================================
# Per-variant sharing core (underlies 3b, 3d, and the flagged file)
# =====================================================================
# Candidate breadth + context list (from the candidate set)
sharing <- cand[, .(
  n_contexts_as_candidate = uniqueN(ctx),
  contexts_list           = paste(sort(unique(ctx)), collapse = ";"),
  dist                    = dist[1L],
  dist_consistent         = uniqueN(dist) == 1L
), by = .(gene, snp)]
if (any(!sharing$dist_consistent))
  log("  WARNING: dist not constant across contexts for",
      sum(!sharing$dist_consistent), "gene x snp (used first).")

# FIMO breadth (from the join)
fimo_ctx <- joined[, .(n_contexts_fimo_hit = uniqueN(ctx)), by = .(gene, snp)]

# Concordance breadth + best primary metric across ALL tested contexts
# (mech_dedup, not just candidate contexts). Computing concordance only on the
# candidate set would make n_contexts_concordant == n_contexts_as_candidate by
# construction (candidates are already concordant-filtered), so we measure it
# over every tested context to keep it informative.
cand_keys <- unique(cand[, .(gene, snp)])
mech_cand_all <- mech_dedup[cand_keys, on = .(gene, snp), nomatch = 0]
conc <- mech_cand_all[, .(
  n_contexts_tested     = uniqueN(ctx),
  n_contexts_concordant = uniqueN(ctx[cbp_eqtl_concordant %in% TRUE]),
  best_cbp_pval         = safe_min(get(PRIMARY))
), by = .(gene, snp)]

sharing <- merge(sharing, fimo_ctx, by = c("gene", "snp"), all.x = TRUE)
sharing[is.na(n_contexts_fimo_hit), n_contexts_fimo_hit := 0L]
sharing <- merge(sharing, conc, by = c("gene", "snp"), all.x = TRUE)

# ---- 3b output ----
b_cols <- c("gene", "snp", "n_contexts_tested", "n_contexts_as_candidate",
            "n_contexts_fimo_hit", "n_contexts_concordant",
            "best_cbp_pval", "contexts_list")
sharing_out <- sharing[, ..b_cols]
setorderv(sharing_out, c("n_contexts_as_candidate", "best_cbp_pval"), order = c(-1, 1))
fwrite(sharing_out, file.path(wd, "variant_sharing_summary.txt"), sep = "\t")
log("3b variant_sharing_summary.txt:", nrow(sharing_out), "unique gene x snp")

# =====================================================================
# 3c -- TF disruption vocabulary (all candidate-FIMO rows, by TF x family x condition)
# =====================================================================
tf_vocab <- joined[, .(
  n_variants_disrupting = uniqueN(snp),
  n_genes_affected      = uniqueN(gene),
  mean_delta_score      = mean(best_abs_delta_score, na.rm = TRUE),
  frac_ALT_gain         = mean(effect_dir == "ALT_gain", na.rm = TRUE)
), by = .(TF, TF_family, condition)]
setorderv(tf_vocab, c("n_variants_disrupting", "mean_delta_score"), order = c(-1, -1))
fwrite(tf_vocab, file.path(wd, "tf_disruption_vocabulary.txt"), sep = "\t")
log("3c tf_disruption_vocabulary.txt:", nrow(tf_vocab), "TF x family x condition rows")

# =====================================================================
# 3d -- TSS distance analysis (all candidate variants, binned)
# =====================================================================
dd <- copy(sharing)
dd[, abs_dist := abs(dist)]
dd[, dist_bin := cut(abs_dist,
                     breaks = c(-Inf, 1000, 10000, 50000, 100000, Inf),
                     labels = c("<1kb", "1-10kb", "10-50kb", "50-100kb", ">100kb"),
                     right = FALSE)]
tss <- dd[, .(
  n_variants    = .N,
  n_shared      = sum(n_contexts_as_candidate > 1),
  frac_shared   = mean(n_contexts_as_candidate > 1),
  mean_cbp_pval = mean(best_cbp_pval, na.rm = TRUE)
), by = dist_bin]
setkey(tss, dist_bin)   # preserves factor (bin) order
fwrite(tss, file.path(wd, "tss_distance_analysis.txt"), sep = "\t")
log("3d tss_distance_analysis.txt:", nrow(tss), "distance bins")

# =====================================================================
# Interpretation layer -- story_candidates_flagged.txt
# =====================================================================
# Pioneer-disruption strength threshold = top quartile of candidate-FIMO disruption
delta_q75 <- quantile(joined$best_priority_abs_delta_log10p, 0.75, na.rm = TRUE)
log("Pioneer strength threshold (top-quartile best_priority_abs_delta_log10p):",
    round(delta_q75, 4))

pio <- joined[, .(pioneer_disruption =
  any(TF_family %in% PIONEER_FAMS & best_priority_abs_delta_log10p >= delta_q75)),
  by = .(gene, snp)]

# Best overall FIMO hit per variant (across its contexts) for a readable headline TF
setorderv(joined, c("gene", "snp",
                    "best_priority_abs_delta_log10p", "best_abs_delta_score", "TF"),
          order = c(1, 1, -1, -1, 1))
top_hit <- joined[, .SD[1L], by = .(gene, snp)][, .(
  gene, snp,
  top_TF = TF, top_TF_family = TF_family,
  top_change_class = strongest_change_class,
  top_delta_log10p = best_priority_abs_delta_log10p,
  top_hit_context  = ctx)]

# Stimulus specificity (candidate in IFNG but not in PBS)
ss <- cand[, .(in_IFNG = any(condition == "IFNG"),
               in_PBS  = any(condition == "PBS")), by = .(gene, snp)]
ss[, stimulus_specific := in_IFNG & !in_PBS]

flagged <- merge(sharing, ss[, .(gene, snp, stimulus_specific)], by = c("gene", "snp"), all.x = TRUE)
flagged <- merge(flagged, pio,     by = c("gene", "snp"), all.x = TRUE)
flagged <- merge(flagged, top_hit, by = c("gene", "snp"), all.x = TRUE)

flagged[is.na(pioneer_disruption), pioneer_disruption := FALSE]
flagged[, shared_regulator     := n_contexts_as_candidate >= 3]
# candidate => concordant (by construction); FIMO hit => TF expressed (is_expressed all TRUE)
flagged[, near_tss_concordant  := abs(dist) < 5000 & n_contexts_fimo_hit > 0]
flagged[, n_flags := shared_regulator + stimulus_specific + pioneer_disruption + near_tss_concordant]

f_cols <- c("gene", "snp", "dist",
            "n_contexts_as_candidate", "n_contexts_fimo_hit", "n_contexts_concordant",
            "contexts_list", "best_cbp_pval",
            "top_TF", "top_TF_family", "top_change_class", "top_delta_log10p", "top_hit_context",
            "shared_regulator", "stimulus_specific", "pioneer_disruption", "near_tss_concordant",
            "n_flags")
flagged_out <- flagged[, ..f_cols]
setorderv(flagged_out, c("n_flags", "best_cbp_pval"), order = c(-1, 1))
fwrite(flagged_out, file.path(wd, "story_candidates_flagged.txt"), sep = "\t")
log("flagged story_candidates_flagged.txt:", nrow(flagged_out), "unique gene x snp")

# ---------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------
cat("\n================ FLAG SUMMARY ================\n")
cat("shared_regulator (>=3 contexts) :", sum(flagged_out$shared_regulator), "\n")
cat("stimulus_specific (IFNG not PBS):", sum(flagged_out$stimulus_specific), "\n")
cat("pioneer_disruption              :", sum(flagged_out$pioneer_disruption), "\n")
cat("near_tss_concordant (<5kb)      :", sum(flagged_out$near_tss_concordant), "\n")
cat("multi-flag variants (n_flags>=2):", sum(flagged_out$n_flags >= 2), "\n")
cat("=============================================\n")
log("Done.")


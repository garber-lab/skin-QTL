#!/usr/bin/env Rscript
# Plot one variant across all ChromBPNet models.
# Usage: Rscript step6_plot_snp_contrib.R <variant_dir> [head=counts] [logo_flank=30] [track_flank=250]
#   e.g. Rscript step6_plot_snp_contrib.R posterchild/rs9604525 counts 30 250
# Output: <variant_dir>/<rsID>.contrib_<head>.pdf
#
# One column per model, rows:
#   1. observed ATAC (pipeline Tn5-shifted bigWig), +/- track_flank
#   2. predicted cut counts, REF vs ALT (fold-averaged), +/- track_flank
#   3. REF contribution logo, +/- logo_flank
#   4. ALT contribution logo, +/- logo_flank
#   5. ALT - REF contribution per position
#   6. hypothetical contributions, all 4 bases, REF context (which base the model prefers)
#   7. ISM heatmap: log2FC for every substitution, variant tile outlined (if ism.tsv exists)
#   8. motif annotation (if motifs/<rsID>.motif_annotation.tsv exists): Fi-NeMo hits for this head,
#      named from MoDISco; red = gained/stronger with ALT, blue = lost/weaker, grey = unchanged;
#      solid = FIMO agrees (same family, same direction), dashed = Fi-NeMo only
# plus a summary table (mean log2FC, ALT/REF, folds > 0, JSD) under the panels

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(patchwork); library(ggseqlogo); library(gridExtra)
})

args        <- commandArgs(trailingOnly = TRUE)
indir       <- args[1]
head_type   <- if (length(args) >= 2) args[2] else "counts"
logo_flank  <- if (length(args) >= 3) as.integer(args[3]) else 30L
track_flank <- if (length(args) >= 4) as.integer(args[4]) else 250L

rsid <- basename(normalizePath(indir))
f    <- function(x) file.path(indir, paste0(rsid, ".", x))
rd   <- function(x, ...) read_tsv(f(x), show_col_types = FALSE, progress = FALSE, ...)

var     <- rd("variant.tsv", col_names = c("chr", "pos", "ref", "alt", "id"), col_types = "ciccc")
contrib <- rd("contrib_window.tsv") %>% filter(shap_type == head_type, abs(pos_rel) <= logo_flank)
pred    <- rd("pred_profiles.tsv")  %>% filter(abs(pos_rel) <= track_flank)
bw      <- rd("bigwig_tracks.tsv")  %>% filter(track == "observed", abs(pos_rel) <= track_flank)
# depth-normalised if the summarizer wrote it, so all models share one y-axis
bw_cpm  <- "value_cpm" %in% names(bw) && all(is.finite(bw$value_cpm))
bw      <- bw %>% mutate(value_plot = if (bw_cpm) value_cpm else value)
bw_ymax <- max(bw$value_plot) * 1.05
scores  <- rd("scores_by_fold.tsv")
has_ism <- file.exists(f("ism.tsv"))
if (has_ism) ism <- rd("ism.tsv") %>% filter(abs(pos_rel) <= logo_flank)
motif_file <- file.path(indir, "motifs", paste0(rsid, ".motif_annotation.tsv"))
has_motif <- file.exists(motif_file)
if (!has_motif) message("NOTE: no motif annotation at ", motif_file, " -> plotting without the motif row")
MAX_MOTIFS <- 4                                   # per model column
MAX_MEMBERS <- 3                                  # FIMO TFs listed per label, then "+N"
shorten_members <- function(x) {                  # full list stays in motif_annotation.tsv
  if (is.na(x) || x == "") return(x)
  v <- strsplit(x, ",")[[1]]
  if (length(v) <= MAX_MEMBERS) paste(v, collapse = ",")
  else paste0(paste(v[seq_len(MAX_MEMBERS)], collapse = ","), "+", length(v) - MAX_MEMBERS)
}
if (has_motif) {
  mot <- read_tsv(motif_file, show_col_types = FALSE, progress = FALSE) %>%
    filter(head == head_type, plot_default, end_rel >= -logo_flank, start_rel <= logo_flank) %>%
    mutate(rank_score = if_else(direction != "none", 1e6 + abs(delta_imp), imp_REF + imp_ALT),
           sym = case_when(direction == "up" ~ "\u25B2 ", direction == "down" ~ "\u25BC ", TRUE ~ ""),
           extra = case_when(!is.na(members_expressed) & members_expressed != "" ~ members_expressed,
                             !is.na(fimo_best_TF) ~ fimo_best_TF, TRUE ~ NA_character_),
           extra = vapply(extra, shorten_members, character(1), USE.NAMES = FALSE),  # cap TF list
           extra = if_else(!is.na(extra) & extra != label_tf, paste0(" \u00b7 ", extra), ""),
           label = paste0(sym, label_tf, extra)) %>%
    group_by(model) %>% slice_max(rank_score, n = MAX_MOTIFS, with_ties = FALSE) %>% ungroup()
  # label placement: left-aligned at the bracket start, or right-aligned at its end if it would run
  # off the panel; then stack brackets + labels into lanes so nothing overlaps
  bp_per_char <- 1.45 * (2 * logo_flank + 1) / 61
  max_chars <- floor((2 * logo_flank + 1) / bp_per_char)     # label can never exceed the panel width
  mot <- mot %>% mutate(
    label = if_else(nchar(label) > max_chars, paste0(substr(label, 1, max_chars - 1), "\u2026"), label),
    lab_w = nchar(label) * bp_per_char,
    right_align = start_rel - 0.5 + lab_w > logo_flank + 0.5,
    lab_x = if_else(right_align, pmin(end_rel + 0.5, logo_flank + 0.5), pmax(start_rel - 0.5, -logo_flank - 0.5)),
    lab_hjust = if_else(right_align, 1, 0),
    lo = pmin(start_rel, if_else(right_align, lab_x - lab_w, lab_x)),
    hi = pmax(end_rel, if_else(right_align, lab_x, lab_x + lab_w)))
  lanes <- function(d) {
    d <- arrange(d, lo); ends <- numeric(0); d$lane <- NA_integer_
    for (i in seq_len(nrow(d))) {
      l <- which(ends < d$lo[i] - 1)[1]
      if (is.na(l)) { ends <- c(ends, d$hi[i]); l <- length(ends) } else ends[l] <- d$hi[i]
      d$lane[i] <- l
    }
    d
  }
  mot <- mot %>% group_by(model) %>% group_modify(~ lanes(.x)) %>% ungroup()
  n_lanes <- max(c(1, mot$lane))
}

# model order: cell type, then PBS before IFN
pref   <- c("KRT_PBS", "KRT_IFN", "FRB_PBS", "FRB_IFN", "MEL_PBS", "MEL_IFN")
models <- c(intersect(pref, unique(contrib$model)), setdiff(unique(contrib$model), pref))

# per-model summary: use scores_summary.tsv from the summarizer if present, else compute it
if (file.exists(f("scores_summary.tsv"))) {
  sc <- rd("scores_summary.tsv")
} else {
  sc <- scores %>% group_by(model) %>%
    summarise(mean_log2fc = mean(logfc), sd_log2fc = sd(logfc),
              folds_pos = sum(logfc > 0), n_folds = n(), mean_jsd = mean(jsd), .groups = "drop") %>%
    mutate(alt_over_ref = 2^mean_log2fc)
}
sc <- sc %>% mutate(model = factor(model, levels = models)) %>% arrange(model) %>%
  mutate(model = as.character(model),
         label = sprintf("%s\nlog2FC(ALT/REF) = %.2f \u00b1 %.2f\n%d/%d folds > 0",
                         model, mean_log2fc, sd_log2fc, folds_pos, n_folds))
lab <- setNames(sc$label, sc$model)

# ---- shared y ranges within each row (so models are comparable) ----
logo_rng <- contrib %>% group_by(model, allele, pos_rel) %>%
  summarise(pos = sum(pmax(proj, 0)), neg = sum(pmin(proj, 0)), .groups = "drop")
logo_ylim <- c(min(logo_rng$neg), max(logo_rng$pos)) * 1.05

delta <- contrib %>% group_by(model, allele, pos_rel) %>%
  summarise(v = sum(proj), .groups = "drop") %>%
  pivot_wider(names_from = allele, values_from = v) %>%
  mutate(d = ALT - REF)
delta_ylim <- range(c(delta$d, 0)) * 1.05
pred_ymax  <- max(pred$pred) * 1.05

hyp_rng <- contrib %>% filter(allele == "REF") %>% group_by(model, pos_rel) %>%
  summarise(pos = sum(pmax(hyp, 0)), neg = sum(pmin(hyp, 0)), .groups = "drop")
hyp_ylim <- c(min(hyp_rng$neg), max(hyp_rng$pos)) * 1.05
if (has_ism) ism_lim <- c(-1, 1) * max(abs(ism$mean_log2fc))

theme_trk <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 8, hjust = 0.5), legend.position = "none")
snp_line  <- geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3)

logo_mat <- function(d, val = "proj") {
  w <- d %>% arrange(pos_rel, nucleotide) %>% select(pos_rel, nucleotide, value = all_of(val)) %>%
    pivot_wider(names_from = pos_rel, values_from = value)
  m <- as.matrix(w[, -1]); rownames(m) <- w$nucleotide
  m[c("A", "C", "G", "T"), , drop = FALSE]
}

brks <- pretty(c(-logo_flank, logo_flank), n = 5)
brks <- brks[abs(brks) <= logo_flank]

plot_logo <- function(mdl, al) {
  m <- logo_mat(filter(contrib, model == mdl, allele == al))
  suppressMessages(ggseqlogo(m, method = "custom", seq_type = "dna") +
    annotate("rect", xmin = logo_flank + 0.5, xmax = logo_flank + 1.5,
             ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "gold3") +
    scale_x_continuous(breaks = brks + logo_flank + 1, labels = brks, expand = c(0, 0)) +
    coord_cartesian(ylim = logo_ylim) +
    labs(x = NULL, y = if (mdl == models[1]) paste0(al, " (", if (al == "REF") var$ref else var$alt, ")\n", head_type, " contrib") else NULL) +
    theme_trk)
}

plot_hyp <- function(mdl) {
  m <- logo_mat(filter(contrib, model == mdl, allele == "REF"), val = "hyp")
  suppressMessages(ggseqlogo(m, method = "custom", seq_type = "dna") +
    annotate("rect", xmin = logo_flank + 0.5, xmax = logo_flank + 1.5,
             ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "gold3") +
    scale_x_continuous(breaks = brks + logo_flank + 1, labels = brks, expand = c(0, 0)) +
    coord_cartesian(ylim = hyp_ylim) +
    labs(x = NULL, y = if (mdl == models[1]) paste0("hypothetical\n", head_type, " contrib\n(REF context)") else NULL) +
    theme_trk)
}

plot_ism <- function(mdl) {
  d <- filter(ism, model == mdl)
  refs <- d %>% distinct(pos_rel, ref_base)
  ggplot(d, aes(pos_rel, factor(mut_base, levels = rev(c("A", "C", "G", "T"))))) +
    geom_tile(aes(fill = mean_log2fc)) +
    geom_tile(data = refs, aes(pos_rel, factor(ref_base, levels = rev(c("A", "C", "G", "T")))),
              fill = "grey85", inherit.aes = FALSE) +
    geom_tile(data = filter(d, is_variant == 1), fill = NA, colour = "black", linewidth = 0.6) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                         limits = ism_lim, name = "ISM\nlog2FC") +
    scale_x_continuous(breaks = brks, expand = c(0, 0)) +
    labs(x = if (has_motif) NULL else "bp from variant", y = if (mdl == models[1]) "ISM substitution\n(grey = REF base)" else NULL) +
    theme_trk +
    theme(legend.position = if (mdl == models[length(models)]) "right" else "none",
          legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 7), legend.text = element_text(size = 6))
}

plot_motif <- function(mdl) {
  d <- filter(mot, model == mdl)
  cols <- c(up = "firebrick", down = "steelblue", none = "grey55")
  p <- ggplot() + snp_line +
    scale_x_continuous(breaks = brks, expand = c(0, 0)) +
    coord_cartesian(xlim = c(-logo_flank - 0.5, logo_flank + 0.5), ylim = c(-n_lanes - 0.4, -0.2)) +
    labs(x = "bp from variant",
         y = if (mdl == models[1]) "motifs (Fi-NeMo)\n\u25B2/\u25BC with ALT\nsolid = FIMO agrees" else NULL) +
    theme_trk + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
  if (!nrow(d)) return(p + annotate("text", x = 0, y = -n_lanes / 2 - 0.2, label = "no motif hits", size = 2.5, colour = "grey50"))
  d <- d %>% mutate(y = -lane, lt = if_else(evidence == "finemo+fimo", "solid", "dashed"))
  p +
    geom_segment(data = d, aes(x = start_rel - 0.5, xend = end_rel + 0.5, y = y, yend = y,
                               colour = direction, linetype = lt), linewidth = 0.7) +
    geom_segment(data = d, aes(x = start_rel - 0.5, xend = start_rel - 0.5, y = y, yend = y + 0.18, colour = direction), linewidth = 0.5) +
    geom_segment(data = d, aes(x = end_rel + 0.5, xend = end_rel + 0.5, y = y, yend = y + 0.18, colour = direction), linewidth = 0.5) +
    geom_text(data = d, aes(x = lab_x, y = y + 0.38, label = label, colour = direction, hjust = lab_hjust),
              size = 2.3, fontface = "bold") +
    scale_colour_manual(values = cols, guide = "none") +
    scale_linetype_identity()
}

panels <- list()
for (mdl in models) {
  p_obs <- ggplot(filter(bw, model == mdl), aes(pos_rel, value_plot)) +
    geom_area(fill = "grey40") + snp_line +
    coord_cartesian(ylim = c(0, bw_ymax)) +
    labs(title = lab[[mdl]], x = NULL,
         y = if (mdl == models[1]) paste0("observed ATAC\n", if (bw_cpm) "(cuts per million)" else "(raw counts)") else NULL) +
    theme_trk

  p_pred <- ggplot(filter(pred, model == mdl), aes(pos_rel, pred, colour = allele)) +
    geom_line(linewidth = 0.4) + snp_line +
    scale_colour_manual(values = c(REF = "grey30", ALT = "firebrick")) +
    coord_cartesian(ylim = c(0, pred_ymax)) +
    labs(x = NULL, y = if (mdl == models[1]) "predicted counts\n(REF grey, ALT red)" else NULL) +
    theme_trk

  p_d <- ggplot(filter(delta, model == mdl), aes(pos_rel, d, fill = d > 0)) +
    geom_col(width = 0.9) + geom_hline(yintercept = 0, linewidth = 0.2) +
    scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "steelblue")) +
    scale_x_continuous(breaks = brks, expand = c(0, 0.5)) +
    coord_cartesian(ylim = delta_ylim) +
    labs(x = NULL, y = if (mdl == models[1]) "ALT \u2212 REF\ncontrib" else NULL) +
    theme_trk

  col <- list(p_obs, p_pred, plot_logo(mdl, "REF"), plot_logo(mdl, "ALT"), p_d, plot_hyp(mdl))
  if (has_ism) col <- c(col, list(plot_ism(mdl)))
  if (has_motif) col <- c(col, list(plot_motif(mdl)))
  panels <- c(panels, col)
}

tab_df <- sc %>% transmute(
  Model = model,
  `Mean log2FC` = sprintf("%.2f \u00b1 %.2f", mean_log2fc, sd_log2fc),
  `ALT/REF` = sprintf("%.2f\u00d7", alt_over_ref),
  `Folds > 0` = sprintf("%d/%d", folds_pos, n_folds),
  `Mean JSD` = sprintf("%.3f", mean_jsd))
if ("ism_rank" %in% names(sc)) {
  tab_df <- tab_df %>% mutate(
    `ISM log2FC` = sprintf("%.2f", sc$ism_log2fc),
    `ISM rank (substitution)` = sprintf("%d / %d", sc$ism_rank, sc$ism_n),
    `ISM rank (position)` = sprintf("%d / %d", sc$ism_pos_rank, sc$ism_n_pos))
}
if (has_motif) {
  top_mot <- mot %>% filter(direction != "none") %>% group_by(model) %>%
    slice_max(rank_score, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(Model = model, `Top motif change` = paste0(sym, label_tf,
                if_else(extra == "", "", extra), " (", change, ")"))
  tab_df <- tab_df %>% left_join(top_mot, by = "Model") %>%
    mutate(`Top motif change` = coalesce(`Top motif change`, "\u2014"))
}
tab <- tableGrob(tab_df, rows = NULL, theme = ttheme_minimal(base_size = 9))

row_h <- c(1, 1, 1.2, 1.2, 0.8, 1.2, if (has_ism) 0.9, if (has_motif) 0.35 + 0.3 * n_lanes)
grid_panels <- wrap_plots(panels, ncol = length(models), byrow = FALSE, heights = row_h)
fig <- (grid_panels / wrap_elements(full = tab)) +
  plot_layout(heights = c(10, 1.6)) +
  plot_annotation(
    title = sprintf("%s  %s:%d  %s>%s", rsid, var$chr, var$pos, var$ref, var$alt),
    subtitle = sprintf("ChromBPNet (bias-corrected, 5-fold average) | %s-head DeepSHAP | logos \u00b1%d bp, tracks \u00b1%d bp",
                       head_type, logo_flank, track_flank),
    theme = theme(plot.title = element_text(face = "bold")))

out <- f(paste0("contrib_", head_type, ".pdf"))
ggsave(out, fig, width = 3.2 * length(models), height = 14 + (if (has_ism) 1.5 else 0) + (if (has_motif) 1.5 else 0), device = cairo_pdf)
cat("Wrote", out, "\n")

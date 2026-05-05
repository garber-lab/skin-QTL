# ─── fimo_visualization_v2.R ─────────────────────────────────────────────────
# Generates a publication-quality SNP motif summary figure from:
#   fimo        : raw FIMO output with TF column (all rows, unfiltered)
#   fimo2       : fimo + TF_family joined from Lambert
#   tf_ranked   : summarized + expression-annotated TF table
#
# Usage: called after your existing pipeline has produced fimo, fimo2, tf_ranked
# Outputs: one PDF per panel + merged PDF to dir/
# ─────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(data.table)

# ─── USER SETTINGS (these come from your pipeline args) ──────────────────────
# this_snp, this_celltype, this_condition, dir should already be set upstream
# fimo, fimo2, tf_ranked should already be in your environment
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
this_snp <- args[2]
this_celltype <- args[3] 
this_condition <- args[4]
out_prefix <- paste0(dir, "/fimo_summary_", this_snp, "_", this_celltype, "_", this_condition)
load(paste0(out_prefix, "objects.RData"))

dir <- "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/rs838146"
this_snp <- "rs838146"
this_celltype <- "FRB"
this_condition <- "TNF"
out_prefix <- paste0(dir, "/fimo_summary_", this_snp, "_", this_celltype, "_", this_condition)
load(paste0(out_prefix, "objects.RData"))

# ─── 1. DERIVE SEQUENCE CONTEXT FROM REFERENCE GENOME SEQUENCE ───────────────
# Read SNP position and sequences directly from files written by bash
snp_pos_in_seq <- as.integer(readLines(paste0(dir, "/snp_position_in_seq.txt")))

ref_full <- toupper(paste(readLines(paste0(dir, "/SNP_slop100_REF.fa"))[-1], collapse = ""))
alt_full <- toupper(paste(readLines(paste0(dir, "/SNP_slop100_ALT.fa"))[-1], collapse = ""))
ref_full <- strsplit(ref_full, "")[[1]]
alt_full <- strsplit(alt_full, "")[[1]]

# Allele identities at SNP position
ref_allele <- ref_full[snp_pos_in_seq]
alt_allele <- alt_full[snp_pos_in_seq]

cat(sprintf("[CHECK] snp_pos_in_seq=%d  ref_allele=%s  alt_allele=%s\n",
            snp_pos_in_seq, ref_allele, alt_allele))

# Positions relative to SNP = 0
positions <- seq_along(ref_full) - snp_pos_in_seq

# Build sequence df for plotting
seq_df <- bind_rows(
  tibble(allele = paste0("REF (", ref_allele, " allele)"),
         pos    = positions,
         base   = ref_full,
         is_snp = positions == 0),
  tibble(allele = paste0("ALT (", alt_allele, " allele)"),
         pos    = positions,
         base   = alt_full,
         is_snp = positions == 0)
) %>%
  mutate(
    allele     = factor(allele, levels = unique(allele)),
    base_color = case_when(
      is_snp & allele == levels(allele)[1] ~ "#2980B9",
      is_snp & allele == levels(allele)[2] ~ "#C0392B",
      base == "G"                           ~ "#27AE60",
      base == "C"                           ~ "#8E44AD",
      TRUE                                  ~ "#95A5A6"
    )
  )

# derive sequence data
snp_insight_df <- fread(paste0(dir,"/SNP_context_insight.txt"))
subtitle_parts <- c(
  paste0(snp_insight_df$ref_allele, "->", snp_insight_df$alt_allele, " ", snp_insight_df$mut_type),
  if (snp_insight_df$gc_direction != "no_change")
    paste0(snp_insight_df$gc_direction, "s local GC content (", snp_insight_df$gc_ref_pct, "% -> ", snp_insight_df$gc_alt_pct, "%)"),
  if (!is.na(snp_insight_df$homopolymer_detail)) snp_insight_df$homopolymer_detail,
  snp_insight_df$cpg_detail,
  snp_insight_df$palindrome_detail
)
auto_subtitle <- paste(Filter(Negate(is.null), subtitle_parts), collapse=" | ")

# ─── 2. PREPARE TF TABLE DATA ────────────────────────────────────────────────
# Pull only interesting_final TFs, in ranked order
plot_tfs <- tf_ranked %>%
  filter(interesting_final) %>%
  arrange(desc(strongest_change_class == "gained_in_ALT" | strongest_change_class == "stronger_in_ALT"),
          desc(best_priority_abs_delta_log10p)) %>%
  mutate(
    effect_dir = case_when(
      strongest_change_class %in% c("gained_in_ALT","stronger_in_ALT") ~ "ALT_gain",
      TRUE ~ "ALT_loss"
    ),
    effect_label = case_when(
      strongest_change_class == "gained_in_ALT"   ~ "gained",
      strongest_change_class == "stronger_in_ALT" ~ "stronger",
      strongest_change_class == "lost_in_ALT"     ~ "lost",
      strongest_change_class == "weaker_in_ALT"   ~ "weaker",
      TRUE ~ strongest_change_class
    ),
    # sign of delta reflects direction: positive = gain in ALT, negative = loss
    delta_label   = sprintf("%+.2f", best_abs_delta_score *
                              ifelse(effect_dir=="ALT_gain", 1, -1)),
    # plain ASCII — "YES" instead of checkmark Unicode which drops in some PDF devices
    DE_label      = ifelse(is_DE,      "YES", "-"),
    induced_label = ifelse(is_induced, "YES", "-"),
    CPM_label     = sprintf("%.1f", meanCPM),
    TF_ordered    = factor(TF, levels = rev(TF))
  )

# ─── 3. DERIVE MOTIF POSITIONS FROM FIMO2 ────────────────────────────────────
# For each interesting TF, find where its best motif hit sits within the 100bp
# sequence window, then express that position relative to the SNP (pos=0).
#
# The SNP position within each motif is found empirically by comparing the REF
# and ALT matched sequences directly — more robust than coordinate math, which
# can be thrown off by strand (FIMO may report reverse-complement sequences).

motif_positions <- fimo2 %>%
  filter(TF %in% plot_tfs$TF) %>%          # restrict to TFs we're actually plotting
  group_by(TF) %>%
  slice_max(priority_abs_delta_log10p,      # one row per TF: pick the hit with the
            n=1, with_ties=FALSE) %>%       #   largest allelic difference to represent it
  ungroup() %>%
  mutate(
    motif_start_rel = start_REF - snp_pos_in_seq,   # bp upstream of SNP (negative = upstream)
    motif_end_rel   = stop_REF  - snp_pos_in_seq,   # bp downstream of SNP (positive = downstream)
    # -log10(qvalue) for each allele; clip to avoid log(0)
    nlq_REF       = round(-log10(pmax(qvalue_REF, 1e-300)), 1),
    nlq_ALT       = round(-log10(pmax(qvalue_ALT, 1e-300)), 1),
    nlq_REF_label = ifelse(is.na(nlq_REF), "-", sprintf("%.1f", nlq_REF)),
    nlq_ALT_label = ifelse(is.na(nlq_ALT), "-", sprintf("%.1f", nlq_ALT)),
    # raw score labels
    score_REF_label = ifelse(is.na(score_REF), "-", sprintf("%.1f", score_REF)),
    score_ALT_label = ifelse(is.na(score_ALT), "-", sprintf("%.1f", score_ALT)),
    # find the differing position directly by comparing REF and ALT strings —
    # more robust than coordinate math, which can be thrown off by strand
    snp_offset1 = mapply(function(r, a) {
      r_chars <- strsplit(r, "")[[1]]
      a_chars <- strsplit(a, "")[[1]]
      which(r_chars != a_chars)[1]   # 1-based index of first mismatch
    }, matched_sequence_REF, matched_sequence_ALT),
    # build REF display: left,SNP_base,right using empirically found position
    motif_display_REF = paste0(
      substr(matched_sequence_REF, 1,               snp_offset1 - 1),
      ",",
      substr(matched_sequence_REF, snp_offset1,     snp_offset1),
      ",",
      substr(matched_sequence_REF, snp_offset1 + 1, nchar(matched_sequence_REF))
    ),
    # build ALT display: same position, ALT sequence
    motif_display_ALT = paste0(
      substr(matched_sequence_ALT, 1,               snp_offset1 - 1),
      ",",
      substr(matched_sequence_ALT, snp_offset1,     snp_offset1),
      ",",
      substr(matched_sequence_ALT, snp_offset1 + 1, nchar(matched_sequence_ALT))
    )
  )

# Verification: print the differing position and both allele displays
cat("[CHECK] Motif SNP position verification:\n")
motif_positions %>%
  select(TF, matched_sequence_REF, matched_sequence_ALT,
         snp_offset1, motif_display_REF, motif_display_ALT) %>%
  print(n = 10)

motif_positions <- motif_positions %>%
  select(TF, motif_start_rel, motif_end_rel,
         matched_sequence_REF, matched_sequence_ALT,
         nlq_REF, nlq_ALT, nlq_REF_label, nlq_ALT_label,
         score_REF_label, score_ALT_label,
         motif_display_REF, motif_display_ALT,
         snp_offset1)

# attach motif coordinates back to the TF table for use in Panel B
# Strip any previously joined columns to prevent .x/.y duplication on re-source
plot_tfs <- plot_tfs %>%
  select(-any_of(c(
    "motif_start_rel", "motif_end_rel",
    "matched_sequence_REF", "matched_sequence_ALT",
    "nlq_REF", "nlq_ALT", "nlq_REF_label", "nlq_ALT_label",
    "score_REF_label", "score_ALT_label",
    "motif_display_REF", "motif_display_ALT",
    "snp_offset1"
  ))) %>%
  left_join(motif_positions, by = "TF")

# ─── 4. FAMILY COLOR PALETTE ─────────────────────────────────────────────────
# Auto-assign colors to the TF families present in plot_tfs
families_present <- unique(plot_tfs$TF_family)

# Curated colors for common families; fallback to auto palette
family_color_map <- c(
  "Rel"              = "#C0392B",
  "NF-kB"            = "#C0392B",
  "E2F"              = "#2980B9",
  "C2H2_ZF"          = "#8E44AD",
  "C2H2_ZF;AT_hook"  = "#9B59B6",
  "Homeodomain"      = "#16A085",
  "bHLH"             = "#D35400",
  "bZIP"             = "#F39C12",
  "Nuclear_receptor"  = "#1ABC9C",
  "Ets"              = "#2ECC71",
  "Forkhead"         = "#3498DB",
  "IRF"              = "#E74C3C",
  "SMAD"             = "#95A5A6",
  "RFX"              = "#5D6D7E",
  "HSF"              = "#CA6F1E",
  "HMG/Sox"          = "#7D3C98",
  "GATA"             = "#117A65",
  "T-box"            = "#884EA0",
  "STAT"             = "#1A5276",
  "KLF"              = "#6C3483",
  "SP"               = "#6C3483"
)

# For families not in map, generate colors
extra_families <- setdiff(families_present, names(family_color_map))
if (length(extra_families) > 0) {
  extra_colors <- scales::hue_pal()(length(extra_families))
  names(extra_colors) <- extra_families
  family_color_map <- c(family_color_map, extra_colors)
}

plot_tfs <- plot_tfs %>%
  mutate(tf_color = family_color_map[TF_family])

# ─── 5. BUILD FAMILY SUMMARY BOXES ───────────────────────────────────────────
gain_tfs <- plot_tfs %>% filter(effect_dir == "ALT_gain")
loss_tfs <- plot_tfs %>% filter(effect_dir == "ALT_loss")

# Group by TF_family within each direction
format_family_group <- function(df) {
  df %>%
    group_by(TF_family) %>%
    summarise(tfs = paste(TF, collapse=", "), .groups="drop") %>%
    mutate(line = paste0(TF_family, ": ", tfs)) %>%
    pull(line) %>%
    paste(collapse="\n")
}

gain_summary <- format_family_group(gain_tfs)
loss_summary <- format_family_group(loss_tfs)

# ─── 6. COLORS & THEME ───────────────────────────────────────────────────────
col_header   <- "#2C3E50"
col_stronger <- "#C0392B"
col_weaker   <- "#2980B9"

base_theme <- theme_minimal(base_family="sans") +
  theme(plot.background = element_rect(fill="white", color=NA))

# ─── PANEL A: SEQUENCE CONTEXT ───────────────────────────────────────────────
# Trim to a window around the SNP for readability (±20 bp)
window <- 20
seq_df_plot <- seq_df %>%
  filter(pos >= -window & pos <= window)

snp_tile_df <- seq_df_plot %>%
  filter(is_snp) %>%
  mutate(tile_fill = ifelse(as.numeric(allele) == 1, "#D6EAF8", "#FADBD8"))

p_seq <- ggplot(seq_df_plot, aes(x=pos, y=as.numeric(allele), label=base)) +
  # highlight SNP position background — fill mapped per allele row
  geom_tile(data  = snp_tile_df,
            aes(x = pos, y = as.numeric(allele), fill = tile_fill),
            width = 0.9, height = 0.8, alpha = 0.5, inherit.aes = FALSE) +
  scale_fill_identity() +
  # bases as text
  geom_text(aes(color=base_color), size=4.5, fontface="bold", family="mono") +
  # SNP marker line
  geom_vline(xintercept=0, linetype="dashed", color="#BDC3C7", linewidth=0.5) +
  scale_color_identity() +
  scale_y_continuous(
    breaks = 1:2,
    labels = levels(seq_df_plot$allele),
    limits = c(0.4, 2.6)
  ) +
  scale_x_continuous(
    breaks = 0,
    labels = paste0("SNP\n(", ref_allele, "->", alt_allele, ")")
  ) +
  labs(
    title    = paste0("Sequence context at ", this_snp),
    subtitle = auto_subtitle
  ) +
  theme_minimal(base_family="sans") +
  theme(
    plot.background  = element_rect(fill="white", color=NA),
    plot.title       = element_text(face="bold", size=12, color=col_header),
    plot.subtitle    = element_text(size=9, color="#7F8C8D", face="italic"),
    axis.title       = element_blank(),
    axis.text.x      = element_text(size=9, color="#E74C3C", face="bold"),
    axis.text.y      = element_text(size=10, face="bold", hjust=1),
    panel.grid       = element_blank()
  )

# ─── PANEL B: FAMILY SUMMARY BOXES ───────────────────────────────────────────
summary_df <- tibble(
  direction   = c("ALT gained / strengthened", "ALT lost / weakened"),
  content     = c(gain_summary, loss_summary),
  y           = c(1.5, 0.5),
  box_color   = c("#FADBD8","#D6EAF8"),
  title_color = c("#C0392B","#2980B9")
)

p_summary <- ggplot(summary_df) +
  geom_rect(aes(xmin=0.03, xmax=0.97, ymin=y-0.44, ymax=y+0.44, fill=box_color),
            color="#BDC3C7", linewidth=0.5) +
  geom_text(aes(x=0.5, y=y+0.27, label=direction, color=title_color),
            size=3.8, fontface="bold", hjust=0.5) +
  geom_text(aes(x=0.5, y=y-0.08, label=content),
            size=2.8, hjust=0.5, color="#2C3E50", lineheight=1.4) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,2), expand=c(0,0)) +
  labs(title="TF family summary") +
  theme_void(base_family="sans") +
  theme(
    plot.title      = element_text(face="bold", size=11, color=col_header,
                                   margin=margin(b=8)),
    plot.background = element_rect(fill="white", color=NA),
    plot.margin     = margin(10,10,10,10)
  )

# ─── PANEL C: RANKED TABLE (paginated, 10 TFs per page) ──────────────────────
# Columns: TF | Family | Effect | score REF | score ALT | Dscore |
#          -log10q REF | -log10q ALT | DE | Induced | CPM | Motif REF | Motif ALT
n_rows_per_page <- 40
col_positions <- c(TF      = 0.3,
                   Family  = 2.0,
                   Effect  = 3.5,
                   sREF    = 4.7,
                   sALT    = 5.5,
                   Delta   = 6.3,
                   qREF    = 7.2,
                   qALT    = 8.1,
                   DE      = 8.9,
                   Induced = 9.7,
                   CPM     = 10.5,
                   MotifR  = 11.3,
                   MotifA  = 13.5)
col_width <- 18.0

table_df <- plot_tfs %>%
  arrange(desc(effect_dir), desc(best_priority_abs_delta_log10p)) %>%
  mutate(row_rank = row_number())


# split into chunks of 40
table_pages <- split(table_df, ceiling(table_df$row_rank / n_rows_per_page))

make_table_page <- function(df_chunk) {
  n_rows <- nrow(df_chunk)
  
  # re-factor TF_ordered within this chunk so y axis runs 1..n_rows
  df_chunk <- df_chunk %>%
    mutate(TF_ordered = factor(TF, levels = rev(TF)))
  
  # conditional divider line between gain/loss groups within this chunk
  n_gain      <- sum(df_chunk$effect_dir == "ALT_gain")
  hline_layer <- if (n_gain > 0 && n_gain < n_rows)
    geom_hline(yintercept = n_gain + 0.5, linetype = "dashed",
               color = "#BDC3C7", linewidth = 0.5) else NULL
  
  ggplot(df_chunk) +
    # row shading — faint, separates gain/loss groups visually
    geom_rect(aes(xmin = 0, xmax = col_width,
                  ymin = as.numeric(TF_ordered) - 0.5,
                  ymax = as.numeric(TF_ordered) + 0.5,
                  fill = effect_dir), alpha = 0.07) +
    # header bar — taller so two-line headers fit without clipping
    annotate("rect", xmin = 0, xmax = col_width,
             ymin = n_rows + 0.5, ymax = n_rows + 2.0, fill = col_header) +
    # TF name
    geom_text(aes(x = col_positions["TF"], y = as.numeric(TF_ordered),
                  label = TF, color = tf_color),
              size = 3.2, fontface = "bold", hjust = 0) +
    # Family
    geom_text(aes(x = col_positions["Family"], y = as.numeric(TF_ordered),
                  label = TF_family),
              size = 2.8, hjust = 0, color = "#555555") +
    # Effect
    geom_text(aes(x        = col_positions["Effect"], y = as.numeric(TF_ordered),
                  label    = effect_label,
                  color    = ifelse(effect_dir == "ALT_gain", col_stronger, col_weaker),
                  fontface = ifelse(effect_dir == "ALT_gain", "bold", "italic")),
              size = 3, hjust = 0) +
    # Score REF
    geom_text(aes(x     = col_positions["sREF"], y = as.numeric(TF_ordered),
                  label = score_REF_label,
                  color = ifelse(score_REF_label == "-", "#BDC3C7", "#555555")),
              size = 2.6, hjust = 0.5) +
    # Score ALT
    geom_text(aes(x     = col_positions["sALT"], y = as.numeric(TF_ordered),
                  label = score_ALT_label,
                  color = ifelse(score_ALT_label == "-", "#BDC3C7", "#555555")),
              size = 2.6, hjust = 0.5) +
    # Dscore — "D" prefix used; Unicode delta dropped by some PDF devices
    geom_text(aes(x     = col_positions["Delta"], y = as.numeric(TF_ordered),
                  label = delta_label,
                  color = ifelse(effect_dir == "ALT_gain", "#C0392B", "#2980B9")),
              size = 2.8, hjust = 0.5, fontface = "bold") +
    # -log10(q) REF — orange if >2 (q<0.01), grey otherwise
    geom_text(aes(x     = col_positions["qREF"], y = as.numeric(TF_ordered),
                  label = nlq_REF_label,
                  color = ifelse(nlq_REF_label == "-", "#BDC3C7",
                                 ifelse(nlq_REF > 2,          "#E67E22", "#555555"))),
              size = 2.6, hjust = 0.5) +
    # -log10(q) ALT
    geom_text(aes(x     = col_positions["qALT"], y = as.numeric(TF_ordered),
                  label = nlq_ALT_label,
                  color = ifelse(nlq_ALT_label == "-", "#BDC3C7",
                                 ifelse(nlq_ALT > 2,          "#E67E22", "#555555"))),
              size = 2.6, hjust = 0.5) +
    # DE — "YES" in green, "-" in grey; plain ASCII, safe in all PDF devices
    geom_text(aes(x     = col_positions["DE"], y = as.numeric(TF_ordered),
                  label = DE_label,
                  color = ifelse(is_DE, "#27AE60", "#BDC3C7")),
              size = 2.8, hjust = 0.5, fontface = "bold") +
    # Induced
    geom_text(aes(x     = col_positions["Induced"], y = as.numeric(TF_ordered),
                  label = induced_label,
                  color = ifelse(is_induced, "#27AE60", "#BDC3C7")),
              size = 2.8, hjust = 0.5, fontface = "bold") +
    # CPM
    geom_text(aes(x = col_positions["CPM"], y = as.numeric(TF_ordered),
                  label = CPM_label),
              size = 2.8, hjust = 0.5, color = "#555555") +
    # Motif REF: ,SNP, delimited, blue, monospace
    geom_text(aes(x = col_positions["MotifR"], y = as.numeric(TF_ordered),
                  label = motif_display_REF),
              size = 2.5, hjust = 0, color = "#2980B9", family = "mono") +
    # Motif ALT: ,SNP, delimited, red, monospace
    geom_text(aes(x = col_positions["MotifA"], y = as.numeric(TF_ordered),
                  label = motif_display_ALT),
              size = 2.5, hjust = 0, color = "#C0392B", family = "mono") +
    # column headers — centered in the taller header band
    annotate("text", x = col_positions["TF"],      y = n_rows + 1.25, label = "TF",
             color = "white", size = 3.2, fontface = "bold", hjust = 0) +
    annotate("text", x = col_positions["Family"],  y = n_rows + 1.25, label = "Family",
             color = "white", size = 3.2, fontface = "bold", hjust = 0) +
    annotate("text", x = col_positions["Effect"],  y = n_rows + 1.25, label = "Effect",
             color = "white", size = 3.2, fontface = "bold", hjust = 0) +
    annotate("text", x = col_positions["sREF"],    y = n_rows + 1.25, label = "score\nREF",
             color = "white", size = 2.8, fontface = "bold", hjust = 0.5) +
    annotate("text", x = col_positions["sALT"],    y = n_rows + 1.25, label = "score\nALT",
             color = "white", size = 2.8, fontface = "bold", hjust = 0.5) +
    annotate("text", x = col_positions["Delta"],   y = n_rows + 1.25, label = "Dscore",
             color = "white", size = 2.8, fontface = "bold", hjust = 0.5) +
    annotate("text", x = col_positions["qREF"],    y = n_rows + 1.25, label = "-log10q\nREF",
             color = "white", size = 2.8, fontface = "bold", hjust = 0.5) +
    annotate("text", x = col_positions["qALT"],    y = n_rows + 1.25, label = "-log10q\nALT",
             color = "white", size = 2.8, fontface = "bold", hjust = 0.5) +
    annotate("text", x = col_positions["DE"],      y = n_rows + 1.25, label = "DE",
             color = "white", size = 3.2, fontface = "bold", hjust = 0.5) +
    annotate("text", x = col_positions["Induced"], y = n_rows + 1.25, label = "induced",
             color = "white", size = 3.2, fontface = "bold", hjust = 0.5) +
    annotate("text", x = col_positions["CPM"],     y = n_rows + 1.25, label = "CPM",
             color = "white", size = 3.2, fontface = "bold", hjust = 0.5) +
    annotate("text", x = col_positions["MotifR"],  y = n_rows + 1.25, label = "Motif REF\n(,X,)",
             color = "white", size = 2.8, fontface = "bold", hjust = 0) +
    annotate("text", x = col_positions["MotifA"],  y = n_rows + 1.25, label = "Motif ALT\n(,X,)",
             color = "white", size = 2.8, fontface = "bold", hjust = 0) +
    hline_layer +
    scale_color_identity() +
    scale_fill_manual(values = c("ALT_gain" = "#C0392B", "ALT_loss" = "#2980B9"), guide = "none") +
    scale_y_continuous(limits = c(0.4, n_rows + 2.1), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, col_width), expand = c(0, 0)) +
    labs(title = paste0("Key motif changes | ", this_celltype, " ", this_condition)) +
    theme_void(base_family = "sans") +
    theme(
      plot.title      = element_text(face = "bold", size = 11, color = col_header,
                                     margin = margin(b = 5)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin     = margin(10, 10, 10, 10)
    )
}

# render all table pages
p_table_pages <- lapply(table_pages, make_table_page)

# ─── TITLE + FOOTNOTE ────────────────────────────────────────────────────────
p_title <- ggplot() +
  annotate("text", x=0.5, y=0.65,
           label=paste0("SNP Motif Analysis: ", this_snp),
           size=7, fontface="bold", color=col_header, hjust=0.5) +
  annotate("text", x=0.5, y=0.25,
           label=paste0("ALT allele (", alt_allele, ") | ",
                        this_celltype, " ", this_condition, " context | ",
                        nrow(plot_tfs), " expressed TFs with allele-sensitive motifs"),
           size=3.8, color="#7F8C8D", hjust=0.5, fontface="italic") +
  theme_void() +
  theme(plot.background=element_rect(fill="white", color=NA))

p_foot <- ggplot() +
  annotate("text", x=0.5, y=0.5, hjust=0.5, vjust=0.5, size=2.6, color="#7F8C8D",
           fontface="italic",
           label=paste0(
             "Dscore = delta FIMO log-odds (+ = gain in ALT, - = loss in ALT)  |  ",
             "-log10q: orange = q<0.01  |  ",
             "DE = differentially expressed (log2FC>1, padj<0.05)  |  ",
             "induced = upregulated under ", this_condition, " vs PBS  |  ",
             "CPM = mean expression in ", this_celltype, " ", this_condition, "  |  ",
             "Motif: ,X, marks SNP position; REF in blue, ALT in red"
           )) +
  theme_void() +
  theme(plot.background=element_rect(fill="white", color=NA))

# ─── SAVE: one panel per page, merged into a single PDF ──────────────────────
fwrite(plot_tfs, paste0(out_prefix,"_table.txt"), sep="\t", quote=F)

panel_files <- c(
   paste0(out_prefix, "_p1_title.pdf"),
   paste0(out_prefix, "_p2_sequence.pdf"),
   paste0(out_prefix, "_p3_summary.pdf")
 )
 ggsave(panel_files[1], p_title,   width = 14, height = 3, units = "in", device = cairo_pdf)
 ggsave(panel_files[2], p_seq,     width = 14, height = 4, units = "in", device = cairo_pdf)
 ggsave(panel_files[3], p_summary, width = 14, height = 4, units = "in", device = cairo_pdf)

# save each table page
table_files <- mapply(function(p, i) {
  f <- paste0(out_prefix, "_p_table_", i, ".pdf")
  ggsave(f, p, width = 14, height = 6, units = "in", device = cairo_pdf)
  f
}, p_table_pages, seq_along(p_table_pages), SIMPLIFY = TRUE)

foot_file <- paste0(out_prefix, "_p_foot.pdf")
ggsave(foot_file, p_foot, width = 14, height = 2, units = "in", device = cairo_pdf)

all_files <- c(panel_files, table_files, foot_file)

while (!is.null(dev.list())) dev.off()

merged_pdf <- paste0(out_prefix, ".pdf")
pdf(merged_pdf, width = 14, height = 10, onefile = TRUE)
print(p_title)
print(p_seq)
print(p_summary)
for (i in seq_along(p_table_pages)) {
  cat("printing table page", i, "\n")
  print(p_table_pages[[i]])
}
print(p_foot)
dev.off()
cat("Multi-page PDF saved:", merged_pdf, "\n")

# PNG of first table page for quick preview
ggsave(paste0(out_prefix, ".png"),
       p_table_pages[[1]], width = 14, height = 6, units = "in", dpi = 150)
cat("PNG preview saved:", paste0(out_prefix, ".png"), "\n")

### clean up intermediate sequence files
file.remove(file.path(dir, c(
  "snp_position_in_seq.txt",
  "SNP_slop100_ALT.fa",
  "SNP_slop100_REF.fa",
  "SNP_slop100.bed",
  "SNP_slop100.fa",
  "SNP.bed",
  "fimo_objects.RData"
)))

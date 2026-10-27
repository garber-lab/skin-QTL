library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript 04_compile_fimo_output_plots.R <dir> <snp> <celltype> <condition>")
}

dir <- args[1]
this_snp <- args[2]
this_celltype <- args[3]
this_condition <- args[4]

out_prefix <- paste0(dir, "/fimo_summary_", this_snp, "_", this_celltype, "_", this_condition)
objects_file <- paste0(out_prefix, "_objects.RData")

cat(sprintf(
  "[FIMO_PLOT_START] snp=%s celltype=%s condition=%s dir=%s\n",
  this_snp, this_celltype, this_condition, dir
))

options(error = function() {
  cat(sprintf(
    "\n[FIMO_PLOT_ERROR] snp=%s celltype=%s condition=%s dir=%s\n",
    this_snp, this_celltype, this_condition, dir
  ), file = stderr())
  traceback(2)
  q(status = 1)
})

if (!file.exists(objects_file)) {
  stop("Missing objects file: ", objects_file)
}
load(objects_file)

if (!exists("variant_info")) {
  variant_info_file <- file.path(dir, "variant_info.tsv")
  variant_info <- if (file.exists(variant_info_file)) {
    fread(variant_info_file)
  } else {
    data.table(variant_id = this_snp, variant_type = "unknown")
  }
}

read_fasta_seq <- function(path) {
  if (!file.exists(path)) stop("Missing FASTA file: ", path)
  seq <- toupper(paste(readLines(path)[-1], collapse = ""))
  strsplit(seq, "")[[1]]
}

variant_pos_in_seq <- as.integer(readLines(file.path(dir, "snp_position_in_seq.txt")))

ref_full <- read_fasta_seq(file.path(dir, "SNP_slop100_REF.fa"))
alt_full <- read_fasta_seq(file.path(dir, "SNP_slop100_ALT.fa"))

variant_type <- if ("variant_type" %in% names(variant_info)) variant_info$variant_type[1] else "unknown"
variant_ref  <- if ("ref" %in% names(variant_info)) variant_info$ref[1] else ref_full[variant_pos_in_seq]
variant_alt  <- if ("alt" %in% names(variant_info)) variant_info$alt[1] else alt_full[variant_pos_in_seq]

ref_start1 <- if ("ref_start1" %in% names(variant_info)) variant_info$ref_start1[1] else variant_pos_in_seq
ref_end1   <- if ("ref_end1" %in% names(variant_info)) variant_info$ref_end1[1] else variant_pos_in_seq
alt_start1 <- if ("alt_start1" %in% names(variant_info)) variant_info$alt_start1[1] else variant_pos_in_seq
alt_end1   <- if ("alt_end1" %in% names(variant_info)) variant_info$alt_end1[1] else variant_pos_in_seq

is_snv_plot <- variant_type == "SNV" &&
  length(ref_full) == length(alt_full) &&
  nchar(variant_ref) == 1 &&
  nchar(variant_alt) == 1

cat(sprintf(
  "[CHECK] variant_pos_in_seq=%d variant_type=%s REF=%s ALT=%s REF_len=%d ALT_len=%d\n",
  variant_pos_in_seq, variant_type, variant_ref, variant_alt,
  length(ref_full), length(alt_full)
))

seq_df <- bind_rows(
  tibble(
    allele = paste0("REF (", variant_ref, ")"),
    pos = seq_along(ref_full) - variant_pos_in_seq,
    base = ref_full,
    is_variant = seq_along(ref_full) >= ref_start1 & seq_along(ref_full) <= ref_end1
  ),
  tibble(
    allele = paste0("ALT (", variant_alt, ")"),
    pos = seq_along(alt_full) - variant_pos_in_seq,
    base = alt_full,
    is_variant = seq_along(alt_full) >= alt_start1 & seq_along(alt_full) <= alt_end1
  )
) %>%
  mutate(
    allele = factor(allele, levels = unique(allele)),
    base_color = case_when(
      is_variant & allele == levels(allele)[1] ~ "#2980B9",
      is_variant & allele == levels(allele)[2] ~ "#C0392B",
      base == "G" ~ "#27AE60",
      base == "C" ~ "#8E44AD",
      TRUE ~ "#95A5A6"
    )
  )

calc_gc_pct <- function(seq) {
  round(100 * sum(seq %in% c("G", "C")) / length(seq), 1)
}

gc_ref_pct <- calc_gc_pct(ref_full)
gc_alt_pct <- calc_gc_pct(alt_full)

gc_direction <- case_when(
  gc_alt_pct > gc_ref_pct ~ "increase",
  gc_alt_pct < gc_ref_pct ~ "decrease",
  TRUE ~ "no_change"
)

mut_type <- case_when(
  variant_type == "SNV" & paste0(variant_ref, variant_alt) %in% c("AG", "GA", "CT", "TC") ~ "transition",
  variant_type == "SNV" ~ "transversion",
  TRUE ~ variant_type
)

i <- variant_pos_in_seq
ref_local <- paste(ref_full[max(1, i - 1):min(length(ref_full), i + 1)], collapse = "")
alt_local <- paste(alt_full[max(1, i - 1):min(length(alt_full), i + 1)], collapse = "")

ref_has_cpg <- str_detect(ref_local, "CG")
alt_has_cpg <- str_detect(alt_local, "CG")

cpg_detail <- case_when(
  ref_has_cpg & !alt_has_cpg ~ "disrupts CpG",
  !ref_has_cpg & alt_has_cpg ~ "creates CpG",
  ref_has_cpg & alt_has_cpg ~ "CpG retained",
  TRUE ~ "no local CpG"
)

get_run_len <- function(seq, idx) {
  base <- seq[idx]
  left <- idx
  right <- idx
  while (left > 1 && seq[left - 1] == base) left <- left - 1
  while (right < length(seq) && seq[right + 1] == base) right <- right + 1
  right - left + 1
}

ref_run <- get_run_len(ref_full, i)
alt_run <- get_run_len(alt_full, i)

homopolymer_detail <- ifelse(
  max(ref_run, alt_run) >= 4,
  paste0("homopolymer run ", ref_run, "->", alt_run),
  NA_character_
)

revcomp <- function(x) {
  chartr("ACGT", "TGCA", paste(rev(strsplit(x, "")[[1]]), collapse = ""))
}

w <- 10
ref_win <- paste(ref_full[max(1, i - w):min(length(ref_full), i + w)], collapse = "")
alt_win <- paste(alt_full[max(1, i - w):min(length(alt_full), i + w)], collapse = "")

palindrome_detail <- case_when(
  ref_win == revcomp(ref_win) & alt_win != revcomp(alt_win) ~ "disrupts local palindrome",
  ref_win != revcomp(ref_win) & alt_win == revcomp(alt_win) ~ "creates local palindrome",
  TRUE ~ "no palindrome change"
)

snp_insight_df <- tibble(
  ref_allele = variant_ref,
  alt_allele = variant_alt,
  mut_type = mut_type,
  gc_direction = gc_direction,
  gc_ref_pct = gc_ref_pct,
  gc_alt_pct = gc_alt_pct,
  homopolymer_detail = homopolymer_detail,
  cpg_detail = cpg_detail,
  palindrome_detail = palindrome_detail
)

fwrite(snp_insight_df, file.path(dir, "SNP_context_insight.txt"), sep = "\t")

subtitle_parts <- c(
  paste0(variant_ref, "->", variant_alt, " ", mut_type),
  if (gc_direction != "no_change") paste0(gc_direction, "s local GC content (", gc_ref_pct, "% -> ", gc_alt_pct, "%)"),
  if (!is.na(homopolymer_detail) && homopolymer_detail != "") homopolymer_detail,
  cpg_detail,
  palindrome_detail
)
auto_subtitle <- paste(Filter(Negate(is.null), subtitle_parts), collapse = " | ")

plot_tfs <- tf_ranked %>%
  filter(interesting_final) %>%
  arrange(
    desc(strongest_change_class %in% c("gained_in_ALT", "stronger_in_ALT")),
    desc(best_priority_abs_delta_log10p)
  ) %>%
  mutate(
    effect_dir = case_when(
      strongest_change_class %in% c("gained_in_ALT", "stronger_in_ALT") ~ "ALT_gain",
      TRUE ~ "ALT_loss"
    ),
    effect_label = case_when(
      strongest_change_class == "gained_in_ALT" ~ "gained",
      strongest_change_class == "stronger_in_ALT" ~ "stronger",
      strongest_change_class == "lost_in_ALT" ~ "lost",
      strongest_change_class == "weaker_in_ALT" ~ "weaker",
      TRUE ~ strongest_change_class
    ),
    delta_label = sprintf("%+.2f", best_abs_delta_score * ifelse(effect_dir == "ALT_gain", 1, -1)),
    DE_label = ifelse(is_DE, "YES", "-"),
    induced_label = ifelse(is_induced, "YES", "-"),
    CPM_label = sprintf("%.1f", meanCPM),
    TF_ordered = factor(TF, levels = rev(TF))
  )

if (nrow(plot_tfs) == 0) {
  stop("No interesting_final TFs to plot for this SNP/celltype/condition.")
}

motif_positions <- fimo2 %>%
  filter(TF %in% plot_tfs$TF) %>%
  group_by(TF) %>%
  slice_max(priority_abs_delta_log10p, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    motif_start_rel = coalesce(start_REF, start_ALT) - variant_pos_in_seq,
    motif_end_rel = coalesce(stop_REF, stop_ALT) - variant_pos_in_seq,
    
    nlq_REF = round(-log10(pmax(qvalue_REF, 1e-300)), 1),
    nlq_ALT = round(-log10(pmax(qvalue_ALT, 1e-300)), 1),
    nlq_REF_label = ifelse(is.na(nlq_REF), "-", sprintf("%.1f", nlq_REF)),
    nlq_ALT_label = ifelse(is.na(nlq_ALT), "-", sprintf("%.1f", nlq_ALT)),
    
    score_REF_label = ifelse(is.na(score_REF), "-", sprintf("%.1f", score_REF)),
    score_ALT_label = ifelse(is.na(score_ALT), "-", sprintf("%.1f", score_ALT)),
    
    snp_offset1 = mapply(function(r, a) {
      if (!is_snv_plot || is.na(r) || is.na(a) || nchar(r) != nchar(a)) return(NA_integer_)
      idx <- which(strsplit(r, "")[[1]] != strsplit(a, "")[[1]])
      if (length(idx) == 0) NA_integer_ else idx[1]
    }, matched_sequence_REF, matched_sequence_ALT),
    
    motif_display_REF = mapply(function(r, idx) {
      if (is.na(r)) return("-")
      if (is.na(idx)) return(r)
      paste0(substr(r, 1, idx - 1), ",", substr(r, idx, idx), ",", substr(r, idx + 1, nchar(r)))
    }, matched_sequence_REF, snp_offset1),
    
    motif_display_ALT = mapply(function(a, idx) {
      if (is.na(a)) return("-")
      if (is.na(idx)) return(a)
      paste0(substr(a, 1, idx - 1), ",", substr(a, idx, idx), ",", substr(a, idx + 1, nchar(a)))
    }, matched_sequence_ALT, snp_offset1)
  ) %>%
  select(
    TF, motif_start_rel, motif_end_rel,
    matched_sequence_REF, matched_sequence_ALT,
    nlq_REF, nlq_ALT, nlq_REF_label, nlq_ALT_label,
    score_REF_label, score_ALT_label,
    motif_display_REF, motif_display_ALT,
    snp_offset1
  )

cat("[CHECK] Motif position/display preview:\n")
motif_positions %>%
  select(TF, matched_sequence_REF, matched_sequence_ALT, snp_offset1, motif_display_REF, motif_display_ALT) %>%
  print(n = 10)

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

family_color_map <- c(
  "Rel" = "#C0392B", "NF-kB" = "#C0392B", "E2F" = "#2980B9",
  "C2H2_ZF" = "#8E44AD", "C2H2_ZF;AT_hook" = "#9B59B6",
  "Homeodomain" = "#16A085", "bHLH" = "#D35400", "bZIP" = "#F39C12",
  "Nuclear_receptor" = "#1ABC9C", "Ets" = "#2ECC71",
  "Forkhead" = "#3498DB", "IRF" = "#E74C3C", "SMAD" = "#95A5A6",
  "RFX" = "#5D6D7E", "HSF" = "#CA6F1E", "HMG/Sox" = "#7D3C98",
  "GATA" = "#117A65", "T-box" = "#884EA0", "STAT" = "#1A5276",
  "KLF" = "#6C3483", "SP" = "#6C3483"
)

families_present <- unique(plot_tfs$TF_family)
extra_families <- setdiff(families_present, names(family_color_map))
if (length(extra_families) > 0) {
  extra_colors <- scales::hue_pal()(length(extra_families))
  names(extra_colors) <- extra_families
  family_color_map <- c(family_color_map, extra_colors)
}

plot_tfs <- plot_tfs %>%
  mutate(tf_color = family_color_map[TF_family])

format_family_group <- function(df) {
  if (nrow(df) == 0) return("None")
  df %>%
    group_by(TF_family) %>%
    summarise(tfs = paste(TF, collapse = ", "), .groups = "drop") %>%
    mutate(line = paste0(TF_family, ": ", tfs)) %>%
    pull(line) %>%
    paste(collapse = "\n")
}

gain_summary <- format_family_group(plot_tfs %>% filter(effect_dir == "ALT_gain"))
loss_summary <- format_family_group(plot_tfs %>% filter(effect_dir == "ALT_loss"))

col_header <- "#2C3E50"
col_stronger <- "#C0392B"
col_weaker <- "#2980B9"

window <- 20
seq_df_plot <- seq_df %>%
  filter(pos >= -window & pos <= window)

variant_tile_df <- seq_df_plot %>%
  filter(is_variant) %>%
  mutate(tile_fill = ifelse(as.numeric(allele) == 1, "#D6EAF8", "#FADBD8"))

p_seq <- ggplot(seq_df_plot, aes(x = pos, y = as.numeric(allele), label = base)) +
  geom_tile(
    data = variant_tile_df,
    aes(x = pos, y = as.numeric(allele), fill = tile_fill),
    width = 0.9, height = 0.8, alpha = 0.5, inherit.aes = FALSE
  ) +
  scale_fill_identity() +
  geom_text(aes(color = base_color), size = 4.5, fontface = "bold", family = "mono") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#BDC3C7", linewidth = 0.5) +
  scale_color_identity() +
  scale_y_continuous(breaks = 1:2, labels = levels(seq_df_plot$allele), limits = c(0.4, 2.6)) +
  scale_x_continuous(breaks = 0, labels = paste0("Variant\n(", variant_ref, "->", variant_alt, ")")) +
  labs(title = paste0("Sequence context at ", this_snp), subtitle = auto_subtitle) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 12, color = col_header),
    plot.subtitle = element_text(size = 9, color = "#7F8C8D", face = "italic"),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 9, color = "#E74C3C", face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold", hjust = 1),
    panel.grid = element_blank()
  )

summary_df <- tibble(
  direction = c("ALT gained / strengthened", "ALT lost / weakened"),
  content = c(gain_summary, loss_summary),
  y = c(1.5, 0.5),
  box_color = c("#FADBD8", "#D6EAF8"),
  title_color = c("#C0392B", "#2980B9")
)

p_summary <- ggplot(summary_df) +
  geom_rect(aes(xmin = 0.03, xmax = 0.97, ymin = y - 0.44, ymax = y + 0.44, fill = box_color),
            color = "#BDC3C7", linewidth = 0.5) +
  geom_text(aes(x = 0.5, y = y + 0.27, label = direction, color = title_color),
            size = 3.8, fontface = "bold", hjust = 0.5) +
  geom_text(aes(x = 0.5, y = y - 0.08, label = content),
            size = 2.8, hjust = 0.5, color = "#2C3E50", lineheight = 1.4) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0)) +
  labs(title = "TF family summary") +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 11, color = col_header, margin = margin(b = 8)),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

n_rows_per_page <- 40
col_positions <- c(
  TF = 0.3, Family = 2.0, Effect = 3.5, sREF = 4.7, sALT = 5.5,
  Delta = 6.3, qREF = 7.2, qALT = 8.1, DE = 8.9, Induced = 9.7,
  CPM = 10.5, MotifR = 11.3, MotifA = 13.5
)
col_width <- 18.0

table_df <- plot_tfs %>%
  arrange(desc(effect_dir), desc(best_priority_abs_delta_log10p)) %>%
  mutate(row_rank = row_number())

table_pages <- split(table_df, ceiling(table_df$row_rank / n_rows_per_page))

make_table_page <- function(df_chunk) {
  n_rows <- nrow(df_chunk)
  
  df_chunk <- df_chunk %>%
    mutate(TF_ordered = factor(TF, levels = rev(TF)))
  
  n_gain <- sum(df_chunk$effect_dir == "ALT_gain")
  hline_layer <- if (n_gain > 0 && n_gain < n_rows) {
    geom_hline(yintercept = n_gain + 0.5, linetype = "dashed", color = "#BDC3C7", linewidth = 0.5)
  } else {
    NULL
  }
  
  ggplot(df_chunk) +
    geom_rect(aes(xmin = 0, xmax = col_width,
                  ymin = as.numeric(TF_ordered) - 0.5,
                  ymax = as.numeric(TF_ordered) + 0.5,
                  fill = effect_dir), alpha = 0.07) +
    annotate("rect", xmin = 0, xmax = col_width, ymin = n_rows + 0.5, ymax = n_rows + 2.0, fill = col_header) +
    geom_text(aes(x = col_positions["TF"], y = as.numeric(TF_ordered), label = TF, color = tf_color),
              size = 3.2, fontface = "bold", hjust = 0) +
    geom_text(aes(x = col_positions["Family"], y = as.numeric(TF_ordered), label = TF_family),
              size = 2.8, hjust = 0, color = "#555555") +
    geom_text(aes(x = col_positions["Effect"], y = as.numeric(TF_ordered), label = effect_label,
                  color = ifelse(effect_dir == "ALT_gain", col_stronger, col_weaker),
                  fontface = ifelse(effect_dir == "ALT_gain", "bold", "italic")),
              size = 3, hjust = 0) +
    geom_text(aes(x = col_positions["sREF"], y = as.numeric(TF_ordered), label = score_REF_label,
                  color = ifelse(score_REF_label == "-", "#BDC3C7", "#555555")),
              size = 2.6, hjust = 0.5) +
    geom_text(aes(x = col_positions["sALT"], y = as.numeric(TF_ordered), label = score_ALT_label,
                  color = ifelse(score_ALT_label == "-", "#BDC3C7", "#555555")),
              size = 2.6, hjust = 0.5) +
    geom_text(aes(x = col_positions["Delta"], y = as.numeric(TF_ordered), label = delta_label,
                  color = ifelse(effect_dir == "ALT_gain", "#C0392B", "#2980B9")),
              size = 2.8, hjust = 0.5, fontface = "bold") +
    geom_text(aes(x = col_positions["qREF"], y = as.numeric(TF_ordered), label = nlq_REF_label,
                  color = ifelse(nlq_REF_label == "-", "#BDC3C7", ifelse(nlq_REF > 2, "#E67E22", "#555555"))),
              size = 2.6, hjust = 0.5) +
    geom_text(aes(x = col_positions["qALT"], y = as.numeric(TF_ordered), label = nlq_ALT_label,
                  color = ifelse(nlq_ALT_label == "-", "#BDC3C7", ifelse(nlq_ALT > 2, "#E67E22", "#555555"))),
              size = 2.6, hjust = 0.5) +
    geom_text(aes(x = col_positions["DE"], y = as.numeric(TF_ordered), label = DE_label,
                  color = ifelse(is_DE, "#27AE60", "#BDC3C7")),
              size = 2.8, hjust = 0.5, fontface = "bold") +
    geom_text(aes(x = col_positions["Induced"], y = as.numeric(TF_ordered), label = induced_label,
                  color = ifelse(is_induced, "#27AE60", "#BDC3C7")),
              size = 2.8, hjust = 0.5, fontface = "bold") +
    geom_text(aes(x = col_positions["CPM"], y = as.numeric(TF_ordered), label = CPM_label),
              size = 2.8, hjust = 0.5, color = "#555555") +
    geom_text(aes(x = col_positions["MotifR"], y = as.numeric(TF_ordered), label = motif_display_REF),
              size = 2.5, hjust = 0, color = "#2980B9", family = "mono") +
    geom_text(aes(x = col_positions["MotifA"], y = as.numeric(TF_ordered), label = motif_display_ALT),
              size = 2.5, hjust = 0, color = "#C0392B", family = "mono") +
    annotate("text", x = col_positions["TF"], y = n_rows + 1.25, label = "TF",
             color = "white", size = 3.2, fontface = "bold", hjust = 0) +
    annotate("text", x = col_positions["Family"], y = n_rows + 1.25, label = "Family",
             color = "white", size = 3.2, fontface = "bold", hjust = 0) +
    annotate("text", x = col_positions["Effect"], y = n_rows + 1.25, label = "Effect",
             color = "white", size = 3.2, fontface = "bold", hjust = 0) +
    annotate("text", x = col_positions["sREF"], y = n_rows + 1.25, label = "score\nREF",
             color = "white", size = 2.8, fontface = "bold", hjust = 0.5) +
    annotate("text", x = col_positions["sALT"], y = n_rows + 1.25, label = "score\nALT",
             color = "white", size = 2.8, fontface = "bold", hjust = 0.5) +
    annotate("text", x = col_positions["Delta"], y = n_rows + 1.25, label = "Dscore",
             color = "white", size = 2.8, fontface = "bold", hjust = 0.5) +
    annotate("text", x = col_positions["qREF"], y = n_rows + 1.25, label = "-log10q\nREF",
             color = "white", size = 2.8, fontface = "bold", hjust = 0.5) +
    annotate("text", x = col_positions["qALT"], y = n_rows + 1.25, label = "-log10q\nALT",
             color = "white", size = 2.8, fontface = "bold", hjust = 0.5) +
    annotate("text", x = col_positions["DE"], y = n_rows + 1.25, label = "DE",
             color = "white", size = 3.2, fontface = "bold", hjust = 0.5) +
    annotate("text", x = col_positions["Induced"], y = n_rows + 1.25, label = "induced",
             color = "white", size = 3.2, fontface = "bold", hjust = 0.5) +
    annotate("text", x = col_positions["CPM"], y = n_rows + 1.25, label = "CPM",
             color = "white", size = 3.2, fontface = "bold", hjust = 0.5) +
    annotate("text", x = col_positions["MotifR"], y = n_rows + 1.25,
             label = ifelse(is_snv_plot, "Motif REF\n(,X,)", "Motif REF\n(overlap)"),
             color = "white", size = 2.8, fontface = "bold", hjust = 0) +
    annotate("text", x = col_positions["MotifA"], y = n_rows + 1.25,
             label = ifelse(is_snv_plot, "Motif ALT\n(,X,)", "Motif ALT\n(overlap)"),
             color = "white", size = 2.8, fontface = "bold", hjust = 0) +
    hline_layer +
    scale_color_identity() +
    scale_fill_manual(values = c("ALT_gain" = "#C0392B", "ALT_loss" = "#2980B9"), guide = "none") +
    scale_y_continuous(limits = c(0.4, n_rows + 2.1), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, col_width), expand = c(0, 0)) +
    labs(title = paste0("Key motif changes | ", this_celltype, " ", this_condition)) +
    theme_void(base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 11, color = col_header, margin = margin(b = 5)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
}

p_table_pages <- lapply(table_pages, make_table_page)

p_title <- ggplot() +
  annotate("text", x = 0.5, y = 0.65,
           label = paste0("Variant Motif Analysis: ", this_snp),
           size = 7, fontface = "bold", color = col_header, hjust = 0.5) +
  annotate("text", x = 0.5, y = 0.25,
           label = paste0("ALT allele (", variant_alt, ") | ",
                          this_celltype, " ", this_condition, " context | ",
                          nrow(plot_tfs), " expressed TFs with allele-sensitive motifs"),
           size = 3.8, color = "#7F8C8D", hjust = 0.5, fontface = "italic") +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA))

motif_note <- if (is_snv_plot) {
  "Motif: ,X, marks SNP position; REF in blue, ALT in red"
} else {
  "Motif: shown sequences are FIMO hits overlapping the variant interval; REF in blue, ALT in red"
}

p_foot <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5, size = 2.6, color = "#7F8C8D",
           fontface = "italic",
           label = paste0(
             "Dscore = delta FIMO log-odds (+ = gain in ALT, - = loss in ALT)  |  ",
             "-log10q: orange = q<0.01  |  ",
             "DE = differentially expressed (log2FC>1, padj<0.05)  |  ",
             "induced = upregulated under ", this_condition, " vs PBS  |  ",
             "CPM = mean expression in ", this_celltype, " ", this_condition, "  |  ",
             motif_note
           )) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA))

fwrite(plot_tfs, paste0(out_prefix, "_table.txt"), sep = "\t", quote = FALSE)

panel_files <- c(
  paste0(out_prefix, "_p1_title.pdf"),
  paste0(out_prefix, "_p2_sequence.pdf"),
  paste0(out_prefix, "_p3_summary.pdf")
)

ggsave(panel_files[1], p_title, width = 14, height = 3, units = "in", device = cairo_pdf)
ggsave(panel_files[2], p_seq, width = 14, height = 4, units = "in", device = cairo_pdf)
ggsave(panel_files[3], p_summary, width = 14, height = 4, units = "in", device = cairo_pdf)

table_files <- mapply(function(p, page_i) {
  f <- paste0(out_prefix, "_p_table_", page_i, ".pdf")
  ggsave(f, p, width = 14, height = 6, units = "in", device = cairo_pdf)
  f
}, p_table_pages, seq_along(p_table_pages), SIMPLIFY = TRUE)

foot_file <- paste0(out_prefix, "_p_foot.pdf")
ggsave(foot_file, p_foot, width = 14, height = 2, units = "in", device = cairo_pdf)

while (!is.null(dev.list())) dev.off()

merged_pdf <- paste0(out_prefix, ".pdf")
pdf(merged_pdf, width = 14, height = 10, onefile = TRUE)
print(p_title)
print(p_seq)
print(p_summary)
for (page_i in seq_along(p_table_pages)) {
  cat("printing table page", page_i, "\n")
  print(p_table_pages[[page_i]])
}
print(p_foot)
dev.off()

cat("Multi-page PDF saved:", merged_pdf, "\n")

if (file.exists(merged_pdf)) {
  intermediate_pdfs <- c(
    panel_files,
    table_files,
    foot_file
  )
  intermediate_pdfs <- unique(intermediate_pdfs[file.exists(intermediate_pdfs)])
  if (length(intermediate_pdfs) > 0) {
    file.remove(intermediate_pdfs)
    cat("Removed intermediate panel PDFs:\n")
    cat(paste0("  ", intermediate_pdfs, collapse = "\n"), "\n")
  }
}
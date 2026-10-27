#!/usr/bin/env Rscript
# Combine Fi-NeMo motif hits (model-based) with FIMO (sequence-based), contributions, ISM and
# RNA expression into one motif annotation table per variant.
#
# Usage: Rscript step6_motif_annotate.R <variant_dir> [cb_dir] [bias] [finemo_width=400] [keep_flank=50]
# Needs (all written earlier in the pipeline, under <variant_dir>):
#   motifs/finemo/<model>/<head>/hits.tsv + modisco_source.txt   (Fi-NeMo; region 0 = REF, 1 = ALT)
#   motifs/fimo/<rsID>.fimo_REF_ALT.tsv, <rsID>.fimo_all_hits.tsv (FIMO, optional)
#   <rsID>.contrib_window.tsv, <rsID>.ism.tsv                     (summarizer)
# Output:
#   motifs/<rsID>.motif_annotation.tsv   one row per Fi-NeMo motif instance (REF/ALT paired)
#
# Logic: arrows come from Fi-NeMo only. FIMO hits are used only when they overlap a Fi-NeMo span
# AND belong to the same TF family (Lambert DBD), to name specific TFs and confirm direction.

suppressPackageStartupMessages({
  # data.table / magrittr / tibble first: the RNA helper scripts need them (same set as
  # 04_run_fimo_analyze.R); dplyr & co. loaded after so their functions win name clashes
  library(data.table); library(magrittr); library(tibble)
  library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr)
})

args   <- commandArgs(trailingOnly = TRUE)
vdir   <- normalizePath(args[1])
CB     <- if (length(args) >= 2) args[2] else "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet"
BIAS   <- if (length(args) >= 3) args[3] else "0.5"
WIDTH  <- if (length(args) >= 4) as.integer(args[4]) else 400L
KEEP   <- if (length(args) >= 5) as.integer(args[5]) else 50L   # keep hits within +/- KEEP of the variant
rsid   <- basename(vdir)
center <- WIDTH %/% 2                       # Fi-NeMo window position of the variant (0-based)

############################ CONFIG ############################
LAMBERT <- Sys.getenv("LAMBERT_TABLE", "/pi/manuel.garber-umw/human/skin/eQTLs/literature/Lambert_2018_human_TF_table.csv")
RNA_DIR <- Sys.getenv("RNA_HELPER_DIR", "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025")
RNA_HELPERS <- file.path(RNA_DIR, c("fetch_expressed_genes.R", "fetch_DE_genes_log2FC1_padj0.05.R",
                                    "fetch_induced_genes_log2FC1_padj0.05.R", "fetch_mean_CPM_of_genes.R"))
CHANGE_RATIO <- 1.25                        # ALT/REF coefficient ratio for stronger (1/x for weaker)
COND_MAP <- c(IFN = "IFNG", PBS = "PBS")    # model condition -> RNA-seq condition name
################################################################

f   <- function(x) file.path(vdir, paste0(rsid, ".", x))
rd  <- function(p, ...) read_tsv(p, show_col_types = FALSE, progress = FALSE, ...)
say <- function(...) { cat("[MOTIF]", ..., "\n"); message("[MOTIF] ", paste(...)) }

# ---------------- TF families (Lambert 2018) ----------------
lambert <- if (file.exists(LAMBERT)) {
  read_csv(LAMBERT, show_col_types = FALSE) %>%
    transmute(TF = toupper(Name), family = str_replace_all(DBD, " ", "_")) %>% distinct(TF, .keep_all = TRUE)
} else { say("WARNING: Lambert table not found, families will be NA"); tibble(TF = character(), family = character()) }
fam_of <- function(tf) lambert$family[match(toupper(str_split_fixed(tf, "\\+", 2)[, 1]), lambert$TF)]

# ---------------- MoDISco pattern names from motifs.html ----------------
clean_tf <- function(x) {                   # FOSL1+JUNB_MA1137.1 -> FOSL1+JUNB; Tp73.mouse_p53l_1 -> TP73
  y <- sub("_.*$", "", x); y <- sub("\\.(mouse|human)$", "", y, ignore.case = TRUE); toupper(y)
}
parse_motifs_html <- function(fn) {
  x <- paste(readLines(fn, warn = FALSE), collapse = "\n")
  rows <- str_split(x, "<tr")[[1]][-1]
  cells <- map(rows, ~ str_trim(gsub("<[^>]+>", "", str_match_all(.x, "(?s)<t[hd][^>]*>(.*?)</t[hd]>")[[1]][, 2])))
  hdr <- cells[[1]]
  body <- keep(cells[-1], ~ length(.x) == length(hdr))
  df <- as_tibble(do.call(rbind, body), .name_repair = "minimal"); names(df) <- hdr
  df %>% transmute(pattern, match0 = na_if(match0, ""), qval0 = suppressWarnings(as.numeric(qval0)),
                   num_seqlets = suppressWarnings(as.integer(num_seqlets))) %>%
    mutate(motif_tf = clean_tf(match0), motif_family = fam_of(motif_tf))
}
motifs_html_for <- function(h5) {           # .../fold0/chrombpnet_model/auxiliary/interpret_subsample/x.h5
  file.path(dirname(dirname(dirname(h5))), "evaluation", "modisco_profile", "motifs.html")
}

# ---------------- expression helpers (optional) ----------------
have_rna <- all(file.exists(RNA_HELPERS))
if (have_rna) {
  ok <- tryCatch({ for (s in RNA_HELPERS) source(s); TRUE }, error = function(e) { say("WARNING: RNA helpers failed:", conditionMessage(e)); FALSE })
  have_rna <- ok
} else say("WARNING: RNA helper scripts not found, expression columns will be NA")
expr_sets <- function(ct, cond) {
  if (!have_rna) return(NULL)
  tryCatch(list(
    expressed = unique(fetch_expressed_genes(ct, cond)),
    induced   = if (cond == "PBS") character(0) else unique(fetch_induced_genes(ct, cond)),
    DE        = if (cond == "PBS") character(0) else unique(fetch_DE_genes(ct, cond))),
    error = function(e) { say("WARNING: expression lookup failed for", ct, cond, ":", conditionMessage(e)); NULL })
}

# ---------------- inputs ----------------
contrib <- if (file.exists(f("contrib_window.tsv"))) rd(f("contrib_window.tsv")) else NULL
ism     <- if (file.exists(f("ism.tsv"))) rd(f("ism.tsv")) else NULL
fimo_dir <- file.path(vdir, "motifs", "fimo")
fimo_var <- if (file.exists(file.path(fimo_dir, paste0(rsid, ".fimo_REF_ALT.tsv"))))
  rd(file.path(fimo_dir, paste0(rsid, ".fimo_REF_ALT.tsv"))) else NULL
fimo_all <- if (file.exists(file.path(fimo_dir, paste0(rsid, ".fimo_all_hits.tsv"))))
  rd(file.path(fimo_dir, paste0(rsid, ".fimo_all_hits.tsv"))) else NULL
if (is.null(fimo_var)) say("WARNING: no FIMO table; FIMO support columns will be NA")

fimo_var <- fimo_var %>% { if (is.null(.)) NULL else mutate(.,
  f_start = coalesce(start_rel_ALT, start_rel_REF), f_end = coalesce(stop_rel_ALT, stop_rel_REF),
  f_family = fam_of(TF), f_best_p = pmin(pvalue_REF, pvalue_ALT, na.rm = TRUE)) }
fimo_all <- fimo_all %>% { if (is.null(.)) NULL else mutate(., f_family = fam_of(TF)) }

fin_root <- file.path(vdir, "motifs", "finemo")
models <- if (dir.exists(fin_root)) list.dirs(fin_root, recursive = FALSE, full.names = FALSE) else character(0)
if (!length(models)) { say("WARNING: no Fi-NeMo results under", fin_root, "- nothing to annotate"); quit(status = 0) }

# ---------------- per model x head ----------------
out <- list()
for (m in models) {
  ct   <- str_split_fixed(m, "_", 2)[, 1]
  cond <- unname(COND_MAP[str_split_fixed(m, "_", 2)[, 2]]); if (is.na(cond)) cond <- str_split_fixed(m, "_", 2)[, 2]
  ex   <- expr_sets(ct, cond)

  for (h in list.dirs(file.path(fin_root, m), recursive = FALSE, full.names = FALSE)) {
    hd <- file.path(fin_root, m, h)
    if (!file.exists(file.path(hd, "hits.tsv"))) next
    hits <- rd(file.path(hd, "hits.tsv"))
    if (!nrow(hits)) next

    # pattern names from this model's own MoDISco report
    names_tbl <- tibble(pattern = character(), match0 = character(), qval0 = numeric(),
                        num_seqlets = integer(), motif_tf = character(), motif_family = character())
    src <- file.path(hd, "modisco_source.txt")
    if (file.exists(src)) {
      html <- motifs_html_for(str_trim(readLines(src, warn = FALSE)[1]))
      if (file.exists(html)) names_tbl <- parse_motifs_html(html) else say("WARNING: no motifs.html at", html)
    }

    d <- hits %>% filter(start - center <= KEEP, end - 1L - center >= -KEEP) %>% transmute(
      allele = if_else(peak_id == 0, "REF", "ALT"), pattern = motif_name, strand,
      s = start - center, e = end - 1L - center,                # pos_rel, inclusive (end is exclusive)
      coef = hit_coefficient, sim = hit_similarity, imp = hit_importance) %>%
      mutate(id = row_number())

    # pair REF and ALT instances of the same pattern/strand with overlapping spans
    pr <- inner_join(filter(d, allele == "REF"), filter(d, allele == "ALT"),
                     by = c("pattern", "strand"), suffix = c("_R", "_A"), relationship = "many-to-many") %>%
      filter(s_R <= e_A, s_A <= e_R) %>%
      group_by(id_A) %>% slice_max(coef_R, n = 1, with_ties = FALSE) %>% ungroup() %>%
      group_by(id_R) %>% slice_max(coef_A, n = 1, with_ties = FALSE) %>% ungroup()
    paired <- pr %>% transmute(pattern, strand, start_rel = pmin(s_R, s_A), end_rel = pmax(e_R, e_A),
                               coef_REF = coef_R, coef_ALT = coef_A, sim_REF = sim_R, sim_ALT = sim_A,
                               imp_REF = imp_R, imp_ALT = imp_A)
    lost   <- d %>% filter(allele == "REF", !id %in% pr$id_R) %>%
      transmute(pattern, strand, start_rel = s, end_rel = e, coef_REF = coef, coef_ALT = 0,
                sim_REF = sim, sim_ALT = NA_real_, imp_REF = imp, imp_ALT = 0)
    gained <- d %>% filter(allele == "ALT", !id %in% pr$id_A) %>%
      transmute(pattern, strand, start_rel = s, end_rel = e, coef_REF = 0, coef_ALT = coef,
                sim_REF = NA_real_, sim_ALT = sim, imp_REF = 0, imp_ALT = imp)

    a <- bind_rows(paired, lost, gained) %>%
      mutate(model = m, head = h,
             change = case_when(coef_REF == 0 ~ "gained", coef_ALT == 0 ~ "lost",
                                coef_ALT / coef_REF >= CHANGE_RATIO ~ "stronger",
                                coef_REF / coef_ALT >= CHANGE_RATIO ~ "weaker", TRUE ~ "unchanged"),
             direction = case_when(change %in% c("gained", "stronger") ~ "up",
                                   change %in% c("lost", "weaker") ~ "down", TRUE ~ "none"),
             overlaps_variant = start_rel <= 0 & end_rel >= 0,
             pattern_sign = if_else(str_starts(pattern, "neg"), "neg", "pos"),
             delta_imp = imp_ALT - imp_REF) %>%
      left_join(names_tbl, by = "pattern")

    # model support: contribution over the span (this head) and strongest ISM effect inside it
    a <- a %>% rowwise() %>% mutate(
      contrib_REF = if (is.null(contrib)) NA_real_ else sum(contrib$proj[contrib$model == m & contrib$shap_type == h &
                      contrib$allele == "REF" & between(contrib$pos_rel, start_rel, end_rel)]),
      contrib_ALT = if (is.null(contrib)) NA_real_ else sum(contrib$proj[contrib$model == m & contrib$shap_type == h &
                      contrib$allele == "ALT" & between(contrib$pos_rel, start_rel, end_rel)]),
      ism_max_abs_log2fc = if (is.null(ism)) NA_real_ else {
        v <- ism$mean_log2fc[ism$model == m & between(ism$pos_rel, start_rel, end_rel)]
        if (length(v)) max(abs(v)) else NA_real_ }) %>% ungroup() %>%
      mutate(contrib_delta = contrib_ALT - contrib_REF)

    # FIMO support: same-family HOCOMOCO hits overlapping the span
    a <- a %>% rowwise() %>% mutate(fimo = list({
      out1 <- tibble(fimo_best_TF = NA_character_, fimo_class = NA_character_,
                     fimo_delta_log10p = NA_real_, fimo_members = NA_character_, fimo_supports = NA)
      if (!is.null(fimo_var) && !is.na(motif_family)) {
        fv <- fimo_var %>% filter(f_family == motif_family, f_start <= end_rel, f_end >= start_rel)
        fa <- if (is.null(fimo_all)) NULL else
          fimo_all %>% filter(f_family == motif_family, start_rel <= .env$end_rel, stop_rel >= .env$start_rel)
        members <- unique(c(fv$TF, fa$TF))
        if (length(members)) {
          best <- fv %>% arrange(desc(priority_abs_delta_log10p), f_best_p) %>% slice_head(n = 1)
          if (!nrow(best) && !is.null(fa) && nrow(fa)) best <- fa %>% arrange(pvalue) %>% slice_head(n = 1) %>%
              transmute(TF, motif_change_class = "not_overlapping_variant", delta_log10p_ALT_minus_REF = NA_real_)
          fcls <- best$motif_change_class
          fdir <- case_when(fcls %in% c("gained_in_ALT", "stronger_in_ALT") ~ "up",
                            fcls %in% c("lost_in_ALT", "weaker_in_ALT") ~ "down", TRUE ~ "none")
          out1 <- tibble(fimo_best_TF = best$TF, fimo_class = fcls,
                         fimo_delta_log10p = best$delta_log10p_ALT_minus_REF,
                         fimo_members = paste(members, collapse = ","),
                         fimo_supports = (direction == "none") || (fdir == direction))
        }
      }
      out1 })) %>% ungroup() %>% unnest(fimo)

    # expression of the FIMO-supported family members in this cell type x condition
    a <- a %>% rowwise() %>% mutate(
      members_expressed = if (is.null(ex) || is.na(fimo_members)) NA_character_ else {
        mm <- str_split(fimo_members, ",")[[1]]; mm <- mm[mm %in% ex$expressed]
        if (length(mm)) paste0(mm, if_else(mm %in% ex$induced, "*", ""), collapse = ",") else "" },
      any_member_expressed = if (is.null(ex) || is.na(fimo_members)) NA else
        any(str_split(fimo_members, ",")[[1]] %in% ex$expressed)) %>% ungroup() %>%
      mutate(celltype = ct, rna_condition = cond)

    out[[paste(m, h)]] <- a
  }
}

ann <- bind_rows(out)
if (!nrow(ann)) { say("no motif hits to annotate"); quit(status = 0) }

ann <- ann %>%
  mutate(label_tf = coalesce(motif_tf, pattern),
         plot_default = pattern_sign == "pos",          # negative patterns: table only
         evidence = case_when(!is.na(fimo_supports) & fimo_supports ~ "finemo+fimo",
                              TRUE ~ "finemo_only")) %>%
  select(model, celltype, rna_condition, head, pattern, label_tf, motif_family, match0, qval0, num_seqlets,
         strand, start_rel, end_rel, overlaps_variant, change, direction,
         coef_REF, coef_ALT, sim_REF, sim_ALT, imp_REF, imp_ALT, delta_imp,
         contrib_REF, contrib_ALT, contrib_delta, ism_max_abs_log2fc,
         fimo_best_TF, fimo_class, fimo_delta_log10p, fimo_supports, fimo_members,
         members_expressed, any_member_expressed, evidence, pattern_sign, plot_default) %>%
  arrange(model, head, desc(direction != "none"), desc(abs(delta_imp)))

dir.create(file.path(vdir, "motifs"), showWarnings = FALSE)
outF <- file.path(vdir, "motifs", paste0(rsid, ".motif_annotation.tsv"))
write_tsv(ann, outF, na = "NA")
say("wrote", outF, "(", nrow(ann), "motif instances )")
print(ann %>% filter(pattern_sign == "pos", direction != "none") %>%
        select(model, head, label_tf, change, start_rel, end_rel, coef_REF, coef_ALT, fimo_best_TF, members_expressed),
      n = 40, width = 200)

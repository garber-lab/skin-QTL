#!/usr/bin/env Rscript
# =============================================================================
# step7_make_plot_gene_snp_pair.R  --  paired CPM plots for a gene:variant pair.
#
# For each celltype x stimulation, a PBS-vs-stim panel: each donor is a point per
# condition, a line connects that donor's PBS and stim CPM, faceted by genotype.
#
# Two differences from the old per-celltype plotter:
#   (1) Subtitles pull QTL stats from the LONG results_QC output (celltype x
#       condition x QTLtype), passed in as one concatenated STATS_FILE.
#   (2) Donors are split into the analysis SUBSET (donor_subset_<ct>.txt) vs the
#       EXCLUDED donors (present genome-wide but not in the subset). Subset donors
#       are drawn dark; excluded donors are drawn much lighter. Facet strips list
#       the subset donor counts (and excluded counts) per genotype x condition.
#
# Inputs are pre-sliced tables (small) so this is memory-light:
#   GENO_FILE  : genotype table (CHROM POS ID REF ALT <donor...>), must contain <snp>
#   CPM_FILE   : CPM table (... final_gene <sample...>), must contain <gene>
#   STATS_FILE : long QTL stats (celltype condition QTLtype snp gene ... beta p se)
#
# Usage:
#   Rscript step7_make_plot_gene_snp_pair.R <ROOT> <gene> <snp> \
#           <GENO_FILE> <CPM_FILE> <STATS_FILE> <OUT_PDF> <STIMS(comma)>
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  cat("Usage: Rscript step7_make_plot_gene_snp_pair.R <ROOT> <gene> <snp> <GENO_FILE> <CPM_FILE> <STATS_FILE> <OUT_PDF> <STIMS>\n")
  quit(save = "no", status = 1)
}
ROOT       <- args[[1]]
g          <- args[[2]]
snp        <- args[[3]]
GENO_FILE  <- args[[4]]
CPM_FILE   <- args[[5]]
STATS_FILE <- args[[6]]
OUT_PDF    <- args[[7]]
stim_order <- strsplit(args[[8]], ",", fixed = TRUE)[[1]]

META_FILE  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
ct_order   <- c("FRB","KRT","MEL")
cond_cols  <- c(PBS="gray35", IFNG="#0072B2", IFNB="#E69F00", TNF="#D55E00")

message("Plotting ", g, " : ", snp, "  (stims: ", paste(stim_order, collapse=","), ")")

## ---------- CPM (long) + metadata ----------
cpm <- fread(CPM_FILE)
gcol <- if ("final_gene" %in% names(cpm)) "final_gene" else names(cpm)[1]
drop_cols <- intersect(c("gene","name"), names(cpm))
if (length(drop_cols)) cpm <- cpm[, !drop_cols, with = FALSE]
cpm <- cpm[get(gcol) == g]
if (nrow(cpm) == 0) stop("Gene ", g, " not found in ", CPM_FILE)
cpm <- unique(cpm, by = gcol)
samp_cols <- setdiff(names(cpm), gcol)
expr_long <- melt(cpm, id.vars = gcol, measure.vars = samp_cols,
                  variable.name = "sample", value.name = "CPM")
expr_long[, sample := as.character(sample)]

meta <- fread(META_FILE)
expr_long <- merge(expr_long, meta[, .(sample, donor, condition, celltype)], by = "sample")
expr_long[, `:=`(condition = as.character(condition), celltype = as.character(celltype))]

## ---------- genotype dosage for the SNP ----------
geno_dt <- fread(GENO_FILE)
fixed_cols <- c("CHROM","POS","ID","REF","ALT")
geno_donors <- setdiff(names(geno_dt), fixed_cols)
grow <- geno_dt[ID == snp]
if (nrow(grow) == 0) stop("SNP ", snp, " not found in ", GENO_FILE)
ref <- grow$REF[1]; alt <- grow$ALT[1]
dos <- as.numeric(unlist(grow[1, ..geno_donors]))
geno_df <- data.table(donor = geno_donors, dosage = dos)
geno_df <- geno_df[!is.na(dosage)]
geno_df[, genotype := factor(fifelse(dosage==0, paste0(ref,ref),
                              fifelse(dosage==1, paste0(ref,alt),
                              fifelse(dosage==2, paste0(alt,alt), NA_character_))),
                             levels = c(paste0(ref,ref), paste0(ref,alt), paste0(alt,alt)))]

df <- merge(expr_long, geno_df, by = "donor")

## ---------- long QTL stats for subtitles ----------
stats <- tryCatch(fread(STATS_FILE), error = function(e) data.table())
get_stat <- function(ct, cond, qtl){
  if (nrow(stats) == 0 || !all(c("snp","gene","celltype","condition","QTLtype","beta","p") %in% names(stats)))
    return(list(beta=NA, p=NA))
  idx <- stats$snp == snp & stats$gene == g & stats$celltype == ct &
         stats$condition == cond & stats$QTLtype == qtl
  r <- stats[which(idx)]
  if (nrow(r) == 0) return(list(beta=NA, p=NA))
  list(beta = suppressWarnings(as.numeric(r$beta[1])), p = suppressWarnings(as.numeric(r$p[1])))
}
fmt_stat <- function(s){
  if (is.na(s$beta) || is.na(s$p)) return("beta=NA, p=NA")
  stars <- dplyr::case_when(s$p<1e-5~"***", s$p<1e-4~"**", s$p<1e-3~"*", TRUE~"")
  paste0("beta=", sprintf("%.2f", s$beta), ", p=", format(s$p, digits=2, scientific=TRUE), stars)
}

## ---------- per-celltype subset membership ----------
subset_of_ct <- function(ct){
  f <- file.path(ROOT, "data", paste0("donor_subset_", ct, ".txt"))
  if (!file.exists(f)) return(character(0))
  s <- unique(trimws(readLines(f))); s[nzchar(s)]
}

## ---------- subpanel: PBS vs one stim, in one celltype ----------
make_subpanel <- function(ct, stim){
  sub <- df[celltype == ct & condition %in% c("PBS", stim)]
  if (nrow(sub) == 0) return(NULL)
  sub[, in_subset := donor %in% subset_of_ct(ct)]
  sub[, condition := factor(condition, levels = c("PBS", stim))]

  # donor counts per genotype x condition, split subset / excluded
  cnt <- sub[, .(n = uniqueN(donor)), by = .(genotype, dosage, condition, in_subset)]
  lab <- sub[, .(dosage = dosage[1]), by = genotype]
  gcount <- function(gt, cond, insub) {
    v <- cnt[genotype==gt & condition==cond & in_subset==insub, n]; if (length(v)) v[1] else 0L
  }
  lab[, strip_label := {
    sapply(seq_len(.N), function(i){
      gt <- genotype[i]
      paste0(gt,
             "\nsubset  PBS=", gcount(gt,"PBS",TRUE),  " ", stim, "=", gcount(gt,stim,TRUE),
             "\nexcl    PBS=", gcount(gt,"PBS",FALSE), " ", stim, "=", gcount(gt,stim,FALSE))
    })
  }]
  setorder(lab, dosage)
  lab[, strip_label := factor(strip_label, levels = strip_label)]
  sub <- merge(sub, lab[, .(genotype, strip_label)], by = "genotype")

  sub_ex <- sub[in_subset == FALSE]
  sub_in <- sub[in_subset == TRUE]

  subtitle <- paste0(
    ct, " - PBS vs ", stim, "\n",
    "reQTL ", stim, ": ", fmt_stat(get_stat(ct, stim, "reQTL")), "\n",
    "eQTL PBS: ",        fmt_stat(get_stat(ct, "PBS", "eQTL")),  "\n",
    "eQTL ", stim, ": ", fmt_stat(get_stat(ct, stim, "eQTL"))
  )

  ggplot(mapping = aes(x = condition, y = CPM, group = donor)) +
    # excluded donors -- light
    geom_line(data = sub_ex, color = "grey80", alpha = 0.20, linewidth = 0.3) +
    geom_point(data = sub_ex, color = "grey75", alpha = 0.35, size = 1.4) +
    # subset donors -- dark
    geom_line(data = sub_in, color = "grey55", alpha = 0.35, linewidth = 0.4) +
    geom_point(data = sub_in, aes(color = condition), alpha = 0.95, size = 2.1) +
    scale_color_manual(values = cond_cols) +
    facet_wrap(~ strip_label, nrow = 1) +
    labs(title = subtitle, y = "CPM", x = "") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 7, lineheight = 0.9),
      strip.text = element_text(size = 7, lineheight = 0.85),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 7),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )
}

## ---------- assemble ----------
panels <- list()
for (ct in ct_order) for (stim in stim_order) {
  p <- make_subpanel(ct, stim)
  if (!is.null(p)) panels[[paste(ct, stim)]] <- p
}
if (length(panels) == 0) stop("No panels to plot for ", g, ":", snp)

combined <- wrap_plots(panels, ncol = length(stim_order)) +
  plot_annotation(
    title = paste0(g, " : ", snp, "  [SNP ", grow$CHROM[1], ":", grow$POS[1], "]  (dark=subset, light=excluded donors)")
  ) & theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

message("Writing ", OUT_PDF)
pdf(OUT_PDF, width = 5 * length(stim_order), height = 4 * length(ct_order))
print(combined)
invisible(dev.off())
message("Done")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(magrittr)
  library(readr)
  library(patchwork)
  library(tibble)
})

# ---------------------- Argument Parsing ----------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("
Usage:
  Rscript step6_make_plot_for_pair_in_window.R <ct> <gene> <snp>

Example:
  Rscript step6_make_plot_for_pair_in_window.R MEL ITGA1 rs2548496 \\
", "\n")
  quit(save = "no", status = 1)
}

ct         <- args[[1]]
g <- args[[2]] # ITGA1
snp <- args[[3]] # rs2548496

# # toy example
# ct         <- "MEL"
# g <- "IFT122"
# snp <- "rs61737314"


# ---------------------- Set-up --------------------- ####
basedir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
letter <- substr(g, 1, 1)
dir <- paste0(basedir, "/plots/", letter, "/", g, "/", snp)
GENO_FILE <- paste0(dir,"/genotype.txt")
STATS_FILE <- paste0(dir, "/modelstats.txt")
#OUT_PDF <- paste0(basedir,"/plots/", letter, "/plot_",g,"_",snp,".pdf")
OUT_PDF <- paste0(basedir,"/plots/temp_output/plot_",g,"_",snp,".pdf")

CPM_FILE   <- paste0(dir, "/cpm.txt")
META_FILE  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"

message(paste0("Starting to plot for ", g, " : ", snp))
# ---------------------- Load Expression ---------------------- ####
message("[1] Loading CPM expression...")
cpm <- fread(CPM_FILE) %>%
  select(-c("gene","name")) %>%
  distinct(final_gene, .keep_all = TRUE) %>%
  column_to_rownames("final_gene") %>%
  rownames_to_column("gene_id") 

expr_long <- as.data.table(cpm) %>%
  melt(id.vars = "gene_id", variable.name = "sample", value.name = "CPM")
rm(cpm); invisible(gc())

# ---------------------- Metadata ---------------------- ####
message("[2] Loading metadata...")
meta <- fread(META_FILE)
expr_long <- expr_long %>%
  inner_join(meta, by = "sample") %>%
  mutate(condition = as.character(condition),
         celltype  = as.character(celltype))

# ---------------------- Genotype Dosage ---------------------- ####
message("[3] Loading genotype dosage...")
geno_dt <- fread(GENO_FILE)
fixed_cols <- c("CHROM","POS","ID","REF","ALT")

geno_samples <- setdiff(colnames(geno_dt), fixed_cols)
expr_donors <- unique(expr_long$donor)
keep <- intersect(geno_samples, expr_donors)
geno_dt <- geno_dt[, c(fixed_cols, keep), with=FALSE]
setkey(geno_dt, ID)

dosage_to_label <- function(d, ref, alt){
  factor(ifelse(d==0, paste0(ref,ref),
                ifelse(d==1, paste0(ref,alt),
                       ifelse(d==2, paste0(alt,alt), NA))),
         levels=c(paste0(ref,ref),paste0(ref,alt),paste0(alt,alt)))
}

get_genotype_df <- function(snp){
  row <- geno_dt[.(snp)]
  if (nrow(row) == 0) return(NULL)
  
  ref <- row$REF[1]
  alt <- row$ALT[1]
  
  # Extract numeric dosage vector correctly
  dos <- unlist(row[, ..keep], use.names = TRUE)
  dos <- as.numeric(dos)   # convert "0","1","2" to numeric
  names(dos) <- keep       # make sure donor names are preserved
  
  tibble(
    donor = names(dos),
    genotype = dosage_to_label(dos, ref, alt),
    dosage = dos
  )
}


# ---------------------- Stats ---------------------- ####
message("[4] Loading eQTL and reQTL modeling stats...")
stats <- fread(STATS_FILE, header=TRUE)

get_eQTL_stats <- function(s, g, ct, stim){
  row <- stats %>% filter(snp == s, gene == g, celltype == ct)
  if(nrow(row) == 0) return(list(beta = NA, p = NA))
  list(
    beta = row[[paste0("eQTL_beta_", stim)]][1],
    p    = row[[paste0("eQTL_p_", stim)]][1],
    se   = row[[paste0("eQTL_se_", stim)]][1]
  )
}

get_reQTL_stats <- function(s, g, ct, stim){
  row <- stats %>% filter(snp == s, gene == g, celltype == ct)
  if(nrow(row) == 0) return(list(beta = NA, p = NA))
  list(
    beta = row[[paste0("reQTL_dbeta_", stim)]][1],
    p    = row[[paste0("reQTL_p_", stim)]][1],
    se   = row[[paste0("reQTL_se_", stim)]][1]
  )
}


fmt_stat <- function(beta,p){
  if(is.na(beta)|is.na(p)) return("p = NA, β = NA")
  stars <- case_when(p<1e-5~"***", p<1e-4~"**", p<1e-3~"*", TRUE~"")
  paste0("p=", format(p, digits=1, scientific=TRUE),
         ", beta=", sprintf("%.1f",beta),", ",stars)
}

# ---------------------- Plot Function ----------------------
# ---------------------- Plot Helper (PBS vs condition panel) ----------------------
make_subpanel <- function(df, ct, stim, this_snp, this_gene){
  
  sub <- df %>% filter(celltype == ct, condition %in% c("PBS", stim))
  if(nrow(sub) == 0) return(NULL)
  
  # Count donors per genotype + condition
  counts <- sub %>%
    group_by(genotype, dosage, condition) %>%
    summarize(n = n_distinct(donor), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0) %>%
    mutate(strip_label = paste0(genotype, "\nPBS=", PBS, " ", stim, "=", !!sym(stim)))
  
  # join counts into sub
  sub <- sub %>%
    left_join(counts, by = c("genotype", "dosage"))
  
  # ---- enforce facet order by numeric dosage ----
  sub <- sub %>%
    arrange(dosage) %>%   # ensures 0 → 1 → 2 before converting to factor
    mutate(
      genotype_label = factor(strip_label, levels = unique(strip_label))
    )
  
  # ensure PBS appears left on x-axis
  sub$condition <- factor(sub$condition, levels = c("PBS", stim))
  
  
  # ---- plotting ----
  st.reQTL <- get_reQTL_stats(this_snp, this_gene, ct, stim)
  st.eQTL.PBS <- get_eQTL_stats(this_snp, this_gene, ct, "PBS")
  st.eQTL.stim <- get_eQTL_stats(this_snp, this_gene, ct, stim)
  
  subtitle <- paste0(ct, " — PBS vs ", stim, "\n",
                     "(reQTL ", stim, ": ", fmt_stat(st.reQTL$beta, st.reQTL$p), ")","\n",
                     "(eQTL PBS", ": ", fmt_stat(st.eQTL.PBS$beta, st.eQTL.PBS$p), ")", "\n",
                     "(eQTL ",stim,": ", fmt_stat(st.eQTL.stim$beta, st.eQTL.stim$p), ")")
  
  
  ggplot(sub, aes(x = condition, y = CPM, group = donor)) +
    geom_line(alpha = 0.25, color = "grey70") +
    geom_point(aes(color = condition), size = 2.2, alpha = 0.9) +
    scale_color_manual(values = c(
      "PBS" = "gray35",
      "IFNG" = "#0072B2",
      "IFNB" = "#E69F00",
      "TNF" = "#D55E00"
    )) +
    facet_wrap(~ genotype_label, scales = "fixed", nrow = 1) +
    labs(title = subtitle, y = "CPM", x = "") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 3),
      strip.text = element_text(size = 8, lineheight = 0.8),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      legend.position = "none",
      panel.grid.major.x = element_blank(),   # already removing vertical major grid
      panel.grid.major.y = element_blank(),   # ← removes horizontal major grid
      panel.grid.minor.y = element_blank()    # ← removes horizontal minor grid
    )
  
}


# ---------------------- Main Plot for a SNP : gene ----------------------
make_plot <- function(this_gene, this_snp){
  
  df <- expr_long %>% filter(gene_id == this_gene)
  geno <- get_genotype_df(this_snp)
  df <- df %>% inner_join(geno, by="donor")
  
  ct_order <- c("FRB","KRT","MEL")
  stim_order <- c("IFNG","IFNB","TNF")
  
  panels <- list()
  for(ct in ct_order){
    for(stim in stim_order){
      p <- make_subpanel(df, ct, stim, this_snp, this_gene)
      if(!is.null(p)) panels[[paste(ct, stim)]] <- p
    }
  }
  
  combined <- wrap_plots(panels, ncol=3)
  
  # Apply global title WITHOUT overriding sub-panel title sizes
  combined <- combined + plot_annotation(
    title = paste0(this_gene, ":", this_snp, " [dist ", stats$dist,"bp] [SNP ", geno_dt$CHROM, ":", geno_dt$POS, "]")
  )
  
  combined & theme(
    plot.title = element_text(size=7, face="bold", hjust=0.5)  # affects ONLY main title now
  )
}

# ---------------------- Output ----------------------
message("[5] Writing PDF: ", OUT_PDF)
pdf(OUT_PDF, width=16, height=12)
print(make_plot(g, snp))
dev.off()
message("Done")



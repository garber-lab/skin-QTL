#!/usr/bin/env Rscript

Sys.setenv(TZ = "America/New_York")
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(magrittr)
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(future.apply)
  library(msigdbr)
  library(DESeq2)
})
dir.out = "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/fewIFNBeGene_investigation/"

# =======================
# BIG eGene table was compiled by 
# /pi/manuel.garber-umw/human/skin/eQTLs/DREG/
# 2025models/script/analysis1_cis_eQTL_ATAC_overlap.R
# =======================

##### ----- 1. IFNB eGenes overlap with IFNG/TNF eGenes? ----
dir.in <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap/"
df <- read_tsv(paste0(dir.in, "/cisQTL_sigSNPs_ATACpeak_exprDEinduced_annotated.txt"))

# for each gene x celltype，which in which cytokine conditions there's sig eQTL
# for each celltype x condition, only keep the gene that's induced in that condition.
cytokine_overlap <- df %>%
  filter(condition %in% c("IFNB","IFNG","TNF")) %>%
  filter(
    (condition == "IFNB" & is_induced_FRB_IFNB == TRUE & celltype == "FRB") |
      (condition == "IFNB" & is_induced_KRT_IFNB == TRUE & celltype == "KRT") |
      (condition == "IFNB" & is_induced_MEL_IFNB == TRUE & celltype == "MEL") |
      (condition == "IFNG" & is_induced_FRB_IFNG == TRUE & celltype == "FRB") |
      (condition == "IFNG" & is_induced_KRT_IFNG == TRUE & celltype == "KRT") |
      (condition == "IFNG" & is_induced_MEL_IFNG == TRUE & celltype == "MEL") |
      (condition == "TNF"  & is_induced_FRB_TNF  == TRUE & celltype == "FRB") |
      (condition == "TNF"  & is_induced_KRT_TNF  == TRUE & celltype == "KRT") |
      (condition == "TNF"  & is_induced_MEL_TNF  == TRUE & celltype == "MEL")
  ) %>%
  group_by(gene, celltype) %>%
  summarise(
    has_IFNB = "IFNB" %in% unique(condition),
    has_IFNG = "IFNG" %in% unique(condition),
    has_TNF  = "TNF"  %in% unique(condition),
    .groups  = "drop"
  )

# out of IFNB eGenes, how many of then are IFNG or TNF eGenes?
ifnb_genes <- cytokine_overlap %>%
  filter(has_IFNB == TRUE) %>%
  mutate(
    also_IFNG = has_IFNG,
    also_TNF  = has_TNF,
    IFNB_only = !has_IFNG & !has_TNF
  )

ifnb_genes %>%
  group_by(celltype) %>%
  summarise(
    n_IFNB_eGenes   = n(),
    n_also_IFNG     = sum(also_IFNG),
    n_also_TNF      = sum(also_TNF),
    n_IFNB_only     = sum(IFNB_only),
    pct_IFNB_only   = round(sum(IFNB_only)/n()*100, 1)
  )

# check each cytokine's specificity
# IFNB specificity
ifnb_summary <- cytokine_overlap %>%
  filter(has_IFNB == TRUE) %>%
  group_by(celltype) %>%
  summarise(
    n_eGenes      = n(),
    n_also_IFNG   = sum(has_IFNG),
    n_also_TNF    = sum(has_TNF),
    n_specific    = sum(!has_IFNG & !has_TNF),
    pct_specific  = round(n_specific/n()*100, 1),
    cytokine      = "IFNB"
  )

# IFNG specificity
ifng_summary <- cytokine_overlap %>%
  filter(has_IFNG == TRUE) %>%
  group_by(celltype) %>%
  summarise(
    n_eGenes      = n(),
    n_also_IFNB   = sum(has_IFNB),
    n_also_TNF    = sum(has_TNF),
    n_specific    = sum(!has_IFNB & !has_TNF),
    pct_specific  = round(n_specific/n()*100, 1),
    cytokine      = "IFNG"
  )

# TNF specificity
tnf_summary <- cytokine_overlap %>%
  filter(has_TNF == TRUE) %>%
  group_by(celltype) %>%
  summarise(
    n_eGenes      = n(),
    n_also_IFNB   = sum(has_IFNB),
    n_also_IFNG   = sum(has_IFNG),
    n_specific    = sum(!has_IFNB & !has_IFNG),
    pct_specific  = round(n_specific/n()*100, 1),
    cytokine      = "TNF"
  )

# combined
bind_rows(
  ifnb_summary %>% select(celltype, cytokine, n_eGenes, n_specific, pct_specific),
  ifng_summary %>% select(celltype, cytokine, n_eGenes, n_specific, pct_specific),
  tnf_summary  %>% select(celltype, cytokine, n_eGenes, n_specific, pct_specific)
) %>% arrange(celltype, cytokine)


##### ----- 2. log2FC info of induced eGenes -----
df <- read_tsv(paste0(dir.in, "/cisQTL_sigSNPs_ATACpeak_exprDEinduced_annotated.txt"))
load("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/all_degs_abslog2FC1_padj0.05_post_outlier_exclusion.RData")
vars <- ls(pattern = "^deg")

# compile all DEGs
deg_all <- bind_rows(
  deg.frb.ifnb %>% mutate(celltype = "FRB", condition = "IFNB"),
  deg.frb.ifng %>% mutate(celltype = "FRB", condition = "IFNG"),
  deg.frb.tnf  %>% mutate(celltype = "FRB", condition = "TNF"),
  deg.krt.ifnb %>% mutate(celltype = "KRT", condition = "IFNB"),
  deg.krt.ifng %>% mutate(celltype = "KRT", condition = "IFNG"),
  deg.krt.tnf  %>% mutate(celltype = "KRT", condition = "TNF"),
  deg.mel.ifnb %>% mutate(celltype = "MEL", condition = "IFNB"),
  deg.mel.ifng %>% mutate(celltype = "MEL", condition = "IFNG"),
  deg.mel.tnf  %>% mutate(celltype = "MEL", condition = "TNF")
)

# eQTL genes per celltype x condition (no PBS)
eqtl_genes <- df %>%
  filter(condition != "PBS") %>%
  distinct(gene, celltype, condition)

# overlap
plot_df <- deg_all %>%
  inner_join(eqtl_genes %>% mutate(has_eQTL = TRUE),
             by = c("gene","celltype","condition")) %>%
  mutate(condition = factor(condition, levels = c("IFNB","IFNG","TNF")))

# plot
p <- ggplot(plot_df, aes(x = condition, y = log2FoldChange, fill = condition)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~celltype) +
  scale_fill_manual(values = c(
    "IFNB" = "#4472C4",
    "IFNG" = "#ED7D31",
    "TNF"  = "#70AD47"
  )) +
  labs(
    title    = "log2FoldChange of induced DEGs by cytokine",
    subtitle = "All induced DEGs; eQTL genes overlaid",
    x        = NULL,
    y        = "log2FoldChange",
    fill     = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    legend.position  = "none"
  )

ggsave(paste0(dir.out, "/IFNB_foldchange_comparison.pdf"), p, width = 10, height = 5)


deg_all %>%
  filter(log2FoldChange > 0) %>%  # 只看induced
  group_by(celltype, condition) %>%
  summarise(n_induced = n_distinct(gene))


##### ----- 3. could it be IFNB eQTLs are lower in QTL-beta compared to other cytokines? -------
beta_df <- df %>%
  filter(condition %in% c("IFNB","IFNG","TNF")) %>%
  distinct(snp, gene, celltype, condition, beta) %>%
  mutate(
    abs_beta  = abs(beta),
    condition = factor(condition, levels = c("IFNB","IFNG","TNF"))
  )
n_egenes <- beta_df %>%
  group_by(celltype, condition) %>%
  summarise(n_eGenes = n_distinct(gene), .groups = "drop") %>%
  mutate(label = paste0("n=", n_eGenes))

p_beta <- ggplot(beta_df, aes(x = condition, y = abs_beta, fill = condition)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_text(data = n_egenes,
            aes(x = condition, y = -0.05, label = label),
            size = 3, inherit.aes = FALSE) +
  facet_wrap(~celltype) +
  scale_fill_manual(values = c(
    "IFNB" = "#4472C4",
    "IFNG" = "#ED7D31",
    "TNF"  = "#70AD47"
  )) +
  scale_y_continuous(limits = c(-0.1, 3)) +
  labs(
    title = "eQTL effect size (|beta|) by cytokine condition",
    x     = NULL,
    y     = "|beta|",
    fill  = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    legend.position  = "none"
  )


ggsave(paste0(dir.out, "/eQTL_beta_by_condition.pdf"), p_beta, width = 10, height = 5)


##### ----- 4. could it be IFNB eQTLs have more selection pressure than other cytokines? ----
# in bash, run this:
# workingdir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap/
# SNPbed=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap//QTL_sig_snps_all_celltype_conditions.bed
# prefix=sig_snps_all_celltype_conditions
# bash /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/annotate_SNP_by_overlapping_ATACseq_peaks.sh ${workingdir} ${SNPbed} ${prefix}
phyloP_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap/phyloP_scores_sig_snps_all_celltype_conditions.txt"
plyloP_table <- fread(phyloP_file, header = F) %>%
  set_colnames(c("snp","size","bp_covered_by_phyloP","phyloP_sum","phyloP_mean_in_snp","phyloP_mean_in_covered")) %>%
  dplyr::select(c(snp, phyloP_mean_in_covered)) %>%
  set_colnames(c("snp", "phyloP"))

sig_snps_big_table_with_phyloP <- left_join(sig_snps_big_table, plyloP_table, by="snp")

# filter: cytokine only, induced genes only
phyloP_plot_df <- sig_snps_big_table_with_phyloP %>%
  filter(condition != "PBS") %>%
  filter(
    (condition == "IFNB" & celltype == "FRB" & is_induced_FRB_IFNB == TRUE) |
      (condition == "IFNB" & celltype == "KRT" & is_induced_KRT_IFNB == TRUE) |
      (condition == "IFNB" & celltype == "MEL" & is_induced_MEL_IFNB == TRUE) |
      (condition == "IFNG" & celltype == "FRB" & is_induced_FRB_IFNG == TRUE) |
      (condition == "IFNG" & celltype == "KRT" & is_induced_KRT_IFNG == TRUE) |
      (condition == "IFNG" & celltype == "MEL" & is_induced_MEL_IFNG == TRUE) |
      (condition == "TNF"  & celltype == "FRB" & is_induced_FRB_TNF  == TRUE) |
      (condition == "TNF"  & celltype == "KRT" & is_induced_KRT_TNF  == TRUE) |
      (condition == "TNF"  & celltype == "MEL" & is_induced_MEL_TNF  == TRUE)
  ) %>%
  filter(!is.na(phyloP)) %>%
  mutate(condition = factor(condition, levels = c("IFNB","IFNG","TNF")),
         celltype  = factor(celltype,  levels = c("FRB","KRT","MEL")))

# n SNPs per celltype x condition for annotation
n_snps <- phyloP_plot_df %>%
  group_by(celltype, condition) %>%
  summarise(n = n_distinct(snp), .groups = "drop") %>%
  mutate(label = paste0("n=", n))

# wilcoxon test IFNB vs IFNG and IFNB vs TNF per celltype
stats_df <- phyloP_plot_df %>%
  group_by(celltype) %>%
  summarise(
    p_IFNB_vs_IFNG = wilcox.test(
      phyloP[condition == "IFNB"],
      phyloP[condition == "IFNG"])$p.value,
    p_IFNB_vs_TNF = wilcox.test(
      phyloP[condition == "IFNB"],
      phyloP[condition == "TNF"])$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("IFNB vs IFNG: p=", formatC(p_IFNB_vs_IFNG, format="e", digits=1),
                   "\nIFNB vs TNF: p=",  formatC(p_IFNB_vs_TNF,  format="e", digits=1))
  )

# plot
p_phyloP <- ggplot(phyloP_plot_df, 
                   aes(x = condition, y = phyloP, fill = condition)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(data = n_snps,
            aes(x = condition, y = -8, label = label),
            size = 3, inherit.aes = FALSE) +
  geom_text(data = stats_df,
            aes(x = 2, y = 7, label = label),
            size = 2.8, fontface = "italic", inherit.aes = FALSE) +
  facet_wrap(~celltype) +
  scale_fill_manual(values = c(
    "IFNB" = "#4472C4",
    "IFNG" = "#ED7D31",
    "TNF"  = "#70AD47"
  )) +
  scale_y_continuous(limits = c(-9, 8)) +
  labs(
    title    = "phyloP20way conservation score of eQTL SNPs by cytokine",
    subtitle = "Only induced eGenes per condition; higher = more conserved",
    x        = NULL,
    y        = "phyloP20way score",
    fill     = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    legend.position  = "none"
  )

ggsave(paste0(dir.out, "/phyloP_by_condition.pdf"), p_phyloP, width = 10, height = 5)

##### ----- IFNB low eGene issue: log2FC violin plot [>>> anchor on IFNG eGenes <<<] ####
### load data ----
# big eGene table with annotations
dir.table <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap/"
df <- read_tsv(paste0(dir.table, "/cisQTL_sigSNPs_ATACpeak_exprDEinduced_annotated.txt"))

# compile all DEGs
load("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/all_degs_abslog2FC1_padj0.05_post_outlier_exclusion.RData")
vars <- ls(pattern = "^deg")
deg_all <- bind_rows(
  deg.frb.ifnb %>% mutate(celltype = "FRB", condition = "IFNB"),
  deg.frb.ifng %>% mutate(celltype = "FRB", condition = "IFNG"),
  deg.frb.tnf  %>% mutate(celltype = "FRB", condition = "TNF"),
  deg.krt.ifnb %>% mutate(celltype = "KRT", condition = "IFNB"),
  deg.krt.ifng %>% mutate(celltype = "KRT", condition = "IFNG"),
  deg.krt.tnf  %>% mutate(celltype = "KRT", condition = "TNF"),
  deg.mel.ifnb %>% mutate(celltype = "MEL", condition = "IFNB"),
  deg.mel.ifng %>% mutate(celltype = "MEL", condition = "IFNG"),
  deg.mel.tnf  %>% mutate(celltype = "MEL", condition = "TNF")
)
rm(deg.frb.ifnb, deg.frb.ifng, deg.frb.tnf, 
   deg.krt.ifnb, deg.krt.ifng, deg.krt.tnf,
   deg.mel.ifnb, deg.mel.ifng, deg.mel.tnf)

# compile all DESeq2 results
load("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/all_DESeq2_results_post_outlier_exclusion.RData")
vars <- ls(pattern = "^res.")
deseq_res_all <- bind_rows(
  as.data.frame(res.frb.ifnb) %>% rownames_to_column("gene") %>% mutate(celltype = "FRB", condition = "IFNB"),
  as.data.frame(res.frb.ifng) %>% rownames_to_column("gene") %>% mutate(celltype = "FRB", condition = "IFNG"),
  as.data.frame(res.frb.tnf)  %>% rownames_to_column("gene") %>% mutate(celltype = "FRB", condition = "TNF"),
  as.data.frame(res.krt.ifnb) %>% rownames_to_column("gene") %>% mutate(celltype = "KRT", condition = "IFNB"),
  as.data.frame(res.krt.ifng) %>% rownames_to_column("gene") %>% mutate(celltype = "KRT", condition = "IFNG"),
  as.data.frame(res.krt.tnf)  %>% rownames_to_column("gene") %>% mutate(celltype = "KRT", condition = "TNF"),
  as.data.frame(res.mel.ifnb) %>% rownames_to_column("gene") %>% mutate(celltype = "MEL", condition = "IFNB"),
  as.data.frame(res.mel.ifng) %>% rownames_to_column("gene") %>% mutate(celltype = "MEL", condition = "IFNG"),
  as.data.frame(res.mel.tnf)  %>% rownames_to_column("gene") %>% mutate(celltype = "MEL", condition = "TNF")
)
rm(res.frb.ifnb, res.frb.ifng, res.frb.tnf, 
   res.krt.ifnb, res.krt.ifng, res.krt.tnf,
   res.mel.ifnb, res.mel.ifng, res.mel.tnf)
rm(vars)

### Q4.1: Are IFNG eGenes also induced by IFNB?   ----
celltypes <- c("FRB", "MEL", "KRT")

gene_table_list <- map(celltypes, function(ct) {
  
  # IFNG eGenes for this celltype
  ifng_egenes <- df %>%
    filter(celltype == ct, condition == "IFNG", QTLtype == "eQTL") %>%
    distinct(gene)
  
  # log2FC from deg_all for IFNG and IFNB
  deg_ct <- deg_all %>%
    filter(celltype == ct, condition %in% c("IFNG", "IFNB")) %>%
    select(gene, condition, log2FoldChange)
  
  ifng_egenes %>%
    left_join(deg_ct %>% filter(condition == "IFNG") %>%
                select(gene, log2FC_IFNG = log2FoldChange),
              by = "gene") %>%
    left_join(deg_ct %>% filter(condition == "IFNB") %>%
                select(gene, log2FC_IFNB = log2FoldChange),
              by = "gene") %>%
    mutate(
      celltype = ct,
      induced_in_IFNG = !is.na(log2FC_IFNG),
      induced_in_IFNB = !is.na(log2FC_IFNB)
    )
}) %>%
  set_names(celltypes)

# Quick summary for Q1
map_dfr(celltypes, function(ct) {
  
  ifng_induced <- deg_all %>%
    filter(celltype == ct, condition == "IFNG", padj < 0.05, log2FoldChange > 1.5) %>%
    distinct(gene)
  
  ifng_eg <- df %>%
    filter(celltype == ct, condition == "IFNG", QTLtype == "eQTL") %>%
    distinct(gene)
  
  ifnb_induced <- deg_all %>%
    filter(celltype == ct, condition == "IFNB", padj < 0.05, log2FoldChange > 1.5) %>%
    distinct(gene)
  
  # IFNG eGenes that are also IFNG induced
  ifng_eg_and_induced <- ifng_eg %>%
    filter(gene %in% ifng_induced$gene)
  
  tibble(
    celltype = ct,
    n_IFNG_induced = nrow(ifng_induced),
    n_IFNG_eGenes = nrow(ifng_eg),
    n_IFNG_induced_eGene = nrow(ifng_eg_and_induced),
    # NEW: of those, how many also induced in IFNB?
    n_IFNG_induced_eGene_IFNBinduced = sum(ifng_eg_and_induced$gene %in% ifnb_induced$gene),
    # and as percentage
    pct_also_IFNB = round(100 * sum(ifng_eg_and_induced$gene %in% ifnb_induced$gene) / nrow(ifng_eg_and_induced), 1)
  )
})

# SETUP: define gene sets and log2FC data ----
cytokine_colors <- c(
  "IFNG_sig"     = "#1E8449",   # dark green
  "IFNB_sig"     = "#C0392B",   # dark red
  "IFNB_not_sig" = "#F1948A",   # pink
  "TNF_sig"      = "#1A5276",   # dark blue
  "TNF_not_sig"  = "#85C1E9"    # light blue
)

make_plot_data <- function(anchor_genes_fn) {
  map_dfr(celltypes, function(ct) {
    
    anchor <- anchor_genes_fn(ct)
    
    ifnb_sig_genes <- deg_all %>%
      filter(celltype == ct, condition == "IFNB",
             padj < 0.05, log2FoldChange > 1.5) %>%
      pull(gene)
    
    tnf_sig_genes <- deg_all %>%
      filter(celltype == ct, condition == "TNF",
             padj < 0.05, log2FoldChange > 1.5) %>%
      pull(gene)
    
    deseq_res_all %>%
      filter(celltype == ct,
             condition %in% c("IFNG", "IFNB", "TNF"),  # added IFNG
             gene %in% anchor) %>%
      select(gene, condition, log2FoldChange) %>%
      mutate(
        celltype = ct,
        sig_group = case_when(
          condition == "IFNG"  ~ "IFNG_sig",            # IFNG always sig (it's the anchor)
          condition == "IFNB" & gene %in% ifnb_sig_genes ~ "IFNB_sig",
          condition == "IFNB" & !gene %in% ifnb_sig_genes ~ "IFNB_not_sig",
          condition == "TNF"  & gene %in% tnf_sig_genes  ~ "TNF_sig",
          condition == "TNF"  & !gene %in% tnf_sig_genes ~ "TNF_not_sig"
        )
      )
  })
}

# ANCHOR 1: All IFNG induced genes ----
anchor_ifng_induced <- function(ct) {
  deg_all %>%
    dplyr::filter(celltype == ct, condition == "IFNG",
                  padj < 0.05, log2FoldChange > 1.5) %>%
    pull(gene)
}

# ANCHOR 2: IFNG induced genes that are also IFNG eGenes ----
anchor_ifng_induced_egene <- function(ct) {
  ifng_induced <- deg_all %>%
    dplyr::filter(celltype == ct, condition == "IFNG",
                  padj < 0.05, log2FoldChange > 1.5) %>%
    pull(gene)
  
  ifng_egenes <- df %>%
    dplyr::filter(celltype == ct, condition == "IFNG", QTLtype == "eQTL") %>%
    pull(gene)
  
  intersect(ifng_induced, ifng_egenes)
}

plot_data_induced      <- make_plot_data(anchor_ifng_induced)
plot_data_induced_egene <- make_plot_data(anchor_ifng_induced_egene)


# PLOT FUNCTIONS ----
make_violin <- function(data, title) {
  
  # calculate n for each group
  n_labels <- data %>%
    ungroup() %>%
    dplyr::count(sig_group) %>%
    mutate(
      sig_group = factor(sig_group, levels = names(cytokine_colors)),
      label = paste0("n=", n)
    )
  
  data %>%
    ungroup() %>%
    mutate(sig_group = factor(sig_group, levels = names(cytokine_colors))) %>%
    ggplot(aes(x = sig_group, y = log2FoldChange, fill = sig_group)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 1.5, linetype = "dotted", color = "gray50") +
    geom_text(data = n_labels,
              aes(x = sig_group, label = label),
              y = max(data$log2FoldChange, na.rm = TRUE) * 1.05,
              inherit.aes = FALSE,
              size = 3.5, color = "gray30") +
    scale_fill_manual(values = cytokine_colors) +
    labs(title = title,
         x = NULL, y = "log2FoldChange", fill = NULL) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1))
}

make_scatter <- function(data, compare_condition, title) {
  
  scatter_data <- data %>%
    dplyr::ungroup() %>%
    dplyr::filter(condition %in% c("IFNG", compare_condition)) %>%
    dplyr::select(gene, celltype, condition, log2FoldChange, sig_group) %>%
    tidyr::pivot_wider(names_from = condition,
                       values_from = c(log2FoldChange, sig_group)) %>%
    dplyr::rename(
      log2FC_IFNG       = log2FoldChange_IFNG,
      log2FC_compare    = !!sym(paste0("log2FoldChange_", compare_condition)),
      sig_group_compare = !!sym(paste0("sig_group_", compare_condition))
    )
  
  # unified axis limits
  axis_min <- min(c(scatter_data$log2FC_IFNG,
                    scatter_data$log2FC_compare), na.rm = TRUE) * 1.05
  axis_max <- max(abs(c(scatter_data$log2FC_IFNG,
                        scatter_data$log2FC_compare)), na.rm = TRUE) * 1.05
  axis_lim <- c(axis_min, axis_max)
  
  scatter_data %>%
    dplyr::mutate(sig_group_compare = factor(sig_group_compare,
                                             levels = names(cytokine_colors))) %>%
    ggplot(aes(x = log2FC_IFNG, y = log2FC_compare,
               color = sig_group_compare)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 1.5, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 1.5, linetype = "dotted", color = "gray50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
                color = "gray30", alpha = 0.5) +  # y=x reference line
    scale_color_manual(values = cytokine_colors, drop = FALSE) +
    scale_x_continuous(limits = axis_lim) +
    scale_y_continuous(limits = axis_lim) +
    labs(
      title = title,
      x = "log2FC in IFNG",
      y = paste0("log2FC in ", compare_condition),
      color = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

save_pdf <- function(plot_data, plot_fn, filename, title_prefix) {
  pdf(filename, width = 8, height = 6)
  for (ct in celltypes) {
    p <- plot_fn(
      plot_data %>% filter(celltype == ct),
      title = paste0(title_prefix, " — ", ct)
    )
    print(p)
  }
  dev.off()
  cat("Saved:", filename, "\n")
}

# Violin plots
save_pdf(plot_data_induced,
         make_violin,
         paste0(dir.out,"/violin_IFNG_induced_genes.pdf"),
         "IFNG induced genes: log2FC in IFNB and TNF")

save_pdf(plot_data_induced_egene,
         make_violin,
         paste0(dir.out,"violin_IFNG_induced_eGenes.pdf"),
         "IFNG induced eGenes: log2FC in IFNB and TNF")

# scatter
save_scatter_pdf <- function(plot_data, filename, title_prefix) {
  pdf(filename, width = 7, height = 6)
  for (ct in celltypes) {
    for (cond in c("IFNB", "TNF")) {
      p <- make_scatter(
        data = plot_data %>% dplyr::filter(celltype == ct),
        compare_condition = cond,
        title = paste0(title_prefix, " — ", ct, " vs ", cond)
      )
      print(p)
    }
  }
  dev.off()
  cat("Saved:", filename, "\n")
}

# Save scatter PDFs
save_scatter_pdf(plot_data_induced,
                 paste0(dir.out, "scatter_IFNG_induced_genes.pdf"),
                 "IFNG induced genes: IFNG vs")

save_scatter_pdf(plot_data_induced_egene,
                 paste0(dir.out, "/scatter_IFNG_induced_eGenes.pdf"),
                 "IFNG induced eGenes: IFNG vs")











### Q4.2: For IFNG induced eGenes that are not induced by IFNB, what is their IFNB eQTL beta and pval?
total_res <- data.frame()
for (ct in c("MEL", "KRT", "FRB")) {
  
  # fetch genes ----
  ifng_induced <- deg_all %>%
    dplyr::filter(
      celltype == ct,
      condition == "IFNG",
      padj < 0.05,
      log2FoldChange > 1.5
    ) %>%
    dplyr::pull(gene)
  
  ifnb_induced <- deg_all %>%
    dplyr::filter(
      celltype == ct,
      condition == "IFNB",
      padj < 0.05,
      log2FoldChange > 1.5
    ) %>%
    dplyr::pull(gene)
  
  ifng_egenes <- df %>%
    dplyr::filter(
      celltype == ct,
      condition == "IFNG",
      QTLtype == "eQTL"
    ) %>%
    dplyr::pull(gene)
  
  ifnb_egenes <- df %>%
    dplyr::filter(
      celltype == ct,
      condition == "IFNB",
      QTLtype == "eQTL"
    ) %>%
    dplyr::pull(gene)
  
  ifng_induced_eGenes <- intersect(ifng_induced, ifng_egenes)
  ifng_induced_eGenes_not_ifnb_induced <- setdiff(ifng_induced_eGenes, ifnb_induced)
  
  message(
    ct,
    ", ifng_eGenes: ", length(ifng_egenes),
    ", ifng_induced_eGenes: ", length(ifng_induced_eGenes),
    ", ifng_induced_eGenes_not_ifnb_induced: ", length(ifng_induced_eGenes_not_ifnb_induced)
  )
  
  # fetch IFNB and IFNG eQTL modeling results ----
  dir.res <- paste0(
    "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",
    ct,
    "/data"
  )
  
  chunk_file <- ifelse(
    ct == "MEL",
    paste0(dir.res, "/gene_chunk_dict_old.txt"),
    paste0(dir.res, "/gene_chunk_dict.txt")
  )
  
  chunkF <- fread(chunk_file) %>%
    dplyr::filter(gene %in% ifng_induced_eGenes)
  
  chunkF$chunk <- sprintf("%03d", chunkF$chunk)
  
  total_res_this_celltype <- data.frame()
  
  for (this_gene in chunkF$gene) {
    
    message(this_gene)
    
    # lead SNP as IFNG eQTL of gene
    this_snp <- df %>%
      dplyr::filter(
        celltype == ct,
        condition == "IFNG",
        QTLtype == "eQTL",
        gene == this_gene
      ) %>%
      dplyr::arrange(p) %>%
      dplyr::slice(1) %>%
      dplyr::pull(snp)
    
    chunk <- chunkF %>%
      dplyr::filter(gene == this_gene) %>%
      dplyr::slice(1) %>%
      dplyr::pull(chunk)
    
    this_res <- fread(paste0(dir.res, "/../results/result_", chunk, ".tsv")) %>%
      dplyr::filter(gene == this_gene, snp == this_snp) %>%
      dplyr::select(
        snp,
        gene,
        eQTL_beta_IFNG,
        eQTL_p_IFNG,
        eQTL_beta_IFNB,
        eQTL_p_IFNB,
        qc_n_per_condition_by_genotype
      ) %>%
      dplyr::mutate(celltype = ct) %>%
      tidyr::separate(
        qc_n_per_condition_by_genotype,
        into = c(
          "g0_PBS",  "g0_IFNG", "g0_IFNB", "g0_TNF",
          "g1_PBS",  "g1_IFNG", "g1_IFNB", "g1_TNF",
          "g2_PBS",  "g2_IFNG", "g2_IFNB", "g2_TNF"
        ),
        sep = "_",
        convert = TRUE
      ) %>%
      dplyr::select(
        -contains("PBS"),
        -contains("TNF")
      ) %>%
      dplyr::left_join(
        df %>%
          dplyr::select(-tag, -condition, -QTLtype, -beta, -p),
        by = c("snp", "gene", "celltype")
      ) 
    
    total_res_this_celltype <- dplyr::bind_rows(
      total_res_this_celltype,
      this_res
    )
  }

  
  # define IFNB eGene status ----
  total_res_this_celltype$is_IFNB_eGene <- ifelse(total_res_this_celltype$gene %in% ifnb_egenes, TRUE, FALSE)
  
  # clean-up and rbind to main table ----
  total_res_this_celltype <- total_res_this_celltype %>%
    dplyr::distinct()
  
  total_res <- dplyr::bind_rows(
    total_res,
    total_res_this_celltype
  )
}

# get meanCPM in IFNB and IFNG ----
CPM_FILE      <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"
cpm_matrix <- fread(CPM_FILE) %>%
  dplyr::select(-name, -gene) %>%
  dplyr::filter(final_gene != "ARMCX5-GPRASP2") %>%
  column_to_rownames("final_gene") 

META_FILE     <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"   # columns: sample, donor, condition, etc
meta_matrix <- fread(META_FILE)

cpm_long <- cpm_matrix |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "sample", values_to = "cpm") |>
  left_join(
    meta_matrix |> select(sample, celltype, condition),
    by = "sample"
  ) |>
  group_by(gene, celltype, condition) |>
  summarise(meanCPM = mean(cpm), .groups = "drop") |>
  pivot_wider(names_from = condition,
              values_from = meanCPM,
              names_prefix = "meanCPM_")

total_res <- total_res %>%
  left_join(cpm_long, by = c("gene", "celltype")) %>%
  mutate(
    logCPM_PBS  = log2(meanCPM_PBS  + 1),
    logCPM_IFNB = log2(meanCPM_IFNB + 1),
    logCPM_IFNG = log2(meanCPM_IFNG + 1)
  )


# write to file
fwrite(total_res, file=paste0(dir.out, "/IFNG_inducedEGenes_eQTL_DESeq2_ATAC_feature_matrix.txt"), sep="\t", quote=F)

# IFNG induced eGenes that are not induced in IFNB: scatterplot ----
total_res <- fread(paste0(dir.out, "/IFNG_inducedEGenes_eQTL_DESeq2_ATAC_feature_matrix.txt"))
library(tidyverse)
library(ggrepel)

plot_df <- total_res |>
  rowwise() |>
  mutate(
    induced_IFNG = get(paste0("is_induced_", celltype, "_IFNG")),
    induced_IFNB = get(paste0("is_induced_", celltype, "_IFNB"))
  ) |>
  ungroup() |>
  filter(induced_IFNG == TRUE) |>
  mutate(
    ifnb_status     = if_else(induced_IFNB, "IFNB induced", "not IFNB induced"),
    p_IFNB_neglog10 = -log10(eQTL_p_IFNB),
    label_gene      = gene,
    induced_IFNB    = as.character(induced_IFNB)
  ) |>
  mutate(
    group4 = case_when(
      induced_IFNB == "FALSE" & !is_IFNB_eGene ~ "not induced, not eGene",
      induced_IFNB == "FALSE" &  is_IFNB_eGene ~ "not induced, IS eGene",
      induced_IFNB == "TRUE"  & !is_IFNB_eGene ~ "induced, not eGene",
      induced_IFNB == "TRUE"  &  is_IFNB_eGene ~ "induced, IS eGene"
    )
  )

make_scatter <- function(df, ct) {
  ggplot(df, aes(x     = eQTL_beta_IFNG,
                 y     = eQTL_beta_IFNB,
                 color = group4,
                 shape = group4,
                 size  = p_IFNB_neglog10,
                 label = label_gene)) +
    
    scale_color_manual(values = c(
      "not induced, not eGene" = "grey",  
      "not induced, IS eGene"  = "pink", 
      "induced, not eGene"     = "darkblue", 
      "induced, IS eGene"      = "red" 
    ), name = NULL) +
    
    scale_shape_manual(values = c(
      "not induced, not eGene" = 1,    # 空心圆
      "not induced, IS eGene"  = 16,    # 实心圆
      "induced, not eGene"     = 1,   # 空心圆
      "induced, IS eGene"      = 16    # 实心圆
    ), name = NULL) +
    
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey75") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey75") +
    geom_abline(slope = 1, intercept = 0,
                linewidth = 0.4, linetype = "dashed", color = "grey55") +
    
    geom_point(alpha = 0.8) +
    
    geom_text_repel(
      na.rm         = TRUE,
      size          = 2.8,
      max.overlaps  = 25,
      box.padding   = 0.4,
      segment.color = "grey60",
      segment.size  = 0.3
    ) +
    
    scale_size_continuous(
      range  = c(1.5, 5),
      name   = "-log10(IFNB p)"
    ) +
    
    labs(
      title  = ct,
      x      = "eQTL beta — IFNG lead SNP in IFNG model",
      y      = "eQTL beta — same SNP in IFNB model",
      color  = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      strip.background = element_blank(),
      legend.position  = "bottom",
      legend.box       = "horizontal"
    )
}

celltypes <- c("FRB", "KRT", "MEL")

pdf(paste0(dir.out, "/IFNG_inducedEGenes_eQTL_scatter_byCondition.pdf"),
    width = 12, height = 8)

for (ct in celltypes) {
  df_ct <- filter(plot_df, celltype == ct)
  print(make_scatter(df_ct, ct))
}

dev.off()

# IFNG induced eGenes: decision tree modeling ----
library(rpart)
library(rpart.plot)

# load and compile all DESeq2 results

# prep feature matrix
total_res <- fread(paste0(dir.out, "/IFNG_inducedEGenes_eQTL_DESeq2_ATAC_feature_matrix.txt"))

plot_df <- total_res |>
  rowwise() |>
  mutate(
    induced_IFNG = get(paste0("is_induced_", celltype, "_IFNG")),
    induced_IFNB = get(paste0("is_induced_", celltype, "_IFNB"))
  ) |>
  ungroup() |>
  filter(induced_IFNG == TRUE) |>
  mutate(
    ifnb_status     = if_else(induced_IFNB, "IFNB induced", "not IFNB induced"),
    p_IFNB_neglog10 = -log10(eQTL_p_IFNB),
    label_gene      = gene,
    induced_IFNB    = as.character(induced_IFNB)
  ) |>
  mutate(
    group4 = case_when(
      induced_IFNB == "FALSE" & !is_IFNB_eGene ~ "not induced, not eGene",
      induced_IFNB == "FALSE" &  is_IFNB_eGene ~ "not induced, IS eGene",
      induced_IFNB == "TRUE"  & !is_IFNB_eGene ~ "induced, not eGene",
      induced_IFNB == "TRUE"  &  is_IFNB_eGene ~ "induced, IS eGene"
    )
  )

tree_df <- plot_df |>
  mutate(
    is_IFNB_eGene = as.factor(is_IFNB_eGene),
    g2_diff = g2_IFNG - g2_IFNB
  ) |>
  left_join(
    deseq_res_all |>
      filter(condition == "IFNB") |>
      select(gene, celltype, log2FC_IFNB = log2FoldChange, baseMean_IFNB = baseMean),
    by = c("gene", "celltype")
  ) |>
  left_join(
    deseq_res_all |>
      filter(condition == "IFNG") |>
      select(gene, celltype, log2FC_IFNG = log2FoldChange),
    by = c("gene", "celltype")
  ) |>
  select(
    is_IFNB_eGene,
    # eQTL effect
    eQTL_beta_IFNG,
    eQTL_beta_IFNB,
    eQTL_p_IFNB,
    # genotype power
    g2_IFNG, g2_IFNB, g2_diff,
    # transcription
    induced_IFNB,
    log2FC_IFNB, log2FC_IFNG,
    logCPM_PBS, logCPM_IFNB, logCPM_IFNG,
    # celltype
    celltype
  ) |>
  mutate(
    induced_IFNB = as.factor(induced_IFNB),
    celltype     = as.factor(celltype),
    # replace NA log2FC with 0 (not DE = no fold change)
    log2FC_IFNB  = replace_na(log2FC_IFNB, 0),
    log2FC_IFNG  = replace_na(log2FC_IFNG, 0)
  ) |>
  drop_na()

# quick check
tree_df |> dplyr::count(is_IFNB_eGene)

# remove unecessary features ----
tree_df_bio <- tree_df |> 
  select(-eQTL_p_IFNB, -eQTL_beta_IFNB)

set.seed(42)
tree_fit3 <- rpart(
  is_IFNB_eGene ~ .,
  data    = tree_df_bio,
  method  = "class",
  control = rpart.control(
    maxdepth  = 4,
    minsplit  = 10,
    minbucket = 5,
    cp        = 0.01
  ),
  parms = list(prior = c(
    "FALSE" = sum(tree_df_bio$is_IFNB_eGene == "FALSE") / nrow(tree_df_bio),
    "TRUE"  = sum(tree_df_bio$is_IFNB_eGene == "TRUE")  / nrow(tree_df_bio)
  ))
)

pdf(paste0(dir.out, "/decision_tree_IFNB_eGene_bio_only.pdf"), width = 14, height = 8)
rpart.plot(
  tree_fit3,
  type        = 4,
  extra       = 104,
  box.palette = list("Purples", "Greens"),
  shadow.col  = "grey80",
  nn          = TRUE,
  main        = "Decision tree (bio features only): what predicts IFNB eGene status?"
)
dev.off()

tree_fit3$variable.importance |>
  as.data.frame() |>
  rownames_to_column("feature") |>
  set_colnames(c("feature", "importance")) |>
  arrange(desc(importance)) |>
  ggplot(aes(x = reorder(feature, importance), y = importance)) +
  geom_col(fill = "#4472C4") +
  coord_flip() +
  labs(title = "Feature importance (bio features only)", x = NULL, y = "Importance") +
  theme_classic(base_size = 12)

ggsave(paste0(dir.out, "/decision_tree_feature_importance_bio_only.pdf"), width = 7, height = 5)

# annotate genes based on decision tree nodes ----
# 1. 获取每个基因的叶节点
tree_df_bio$leaf_node <- as.integer(
  rownames(tree_fit3$frame)[tree_fit3$where]
)

# 2. 手动定义 FALSE 叶节点的失败原因
fail_reason_map <- c(
  "2"  = "low genotype power (g2_IFNG < 12)",
  "12" = "low IFNB induction (log2FC_IFNB < 0.77)",
  "28" = "low IFNG expression (logCPM_IFNG < 6.59)"
)

# 3. 标注每个基因
tree_df_annotated <- tree_df_bio |>
  mutate(
    gene     = plot_df$gene,
    fail_reason = case_when(
      is_IFNB_eGene == TRUE             ~ "IS IFNB eGene",
      leaf_node == 2                    ~ fail_reason_map["2"],
      leaf_node == 12                   ~ fail_reason_map["12"],
      leaf_node == 28                   ~ fail_reason_map["28"],
      TRUE                              ~ "other/uncertain"
    )
  )

# 4. C: summary table
tree_df_annotated |>
  filter(is_IFNB_eGene == FALSE) |>
  dplyr::count(fail_reason, celltype) |>
  arrange(fail_reason, celltype) |>
  print()

# bar chart
tree_df_annotated |>
  filter(is_IFNB_eGene == FALSE) |>
  dplyr::count(fail_reason) |>
  ggplot(aes(x = reorder(fail_reason, n), y = n, fill = fail_reason)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(
    "low genotype power (g2_IFNG < 12)"       = "#7F77DD",
    "low IFNB induction (log2FC_IFNB < 0.77)" = "#F4A460",
    "low IFNG expression (logCPM_IFNG < 6.59)"= "#85C1E9",
    "other/uncertain"                          = "grey70"
  )) +
  labs(title = "Why did IFNG eGenes fail to be IFNB eGenes?",
       x = NULL, y = "n genes") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

ggsave(paste0(dir.out, "/fail_reason_barplot.pdf"), width = 8, height = 4)

fwrite(tree_df_annotated, paste0(dir.out,"/IFNG_induced_eGenes_that_arenot_IFNB_eGene_annotated_reason.txt"), sep="\t", quote=F)

# write gene snp pairs for plotting ----
temp <- left_join(tree_df_annotated, total_res[,c("snp","gene","celltype")])

temp1 <- temp[which(temp$celltype=="MEL"),] %>% dplyr::select(gene, snp)
fwrite(temp1, paste0(dir.out,"/temp_MEL_gene_snp_pairs_for_plotting.txt"), quote=F, sep="\t", col.names = F)
temp1 <- temp[which(temp$celltype=="FRB"),] %>% dplyr::select(gene, snp)
fwrite(temp1, paste0(dir.out,"/temp_FRB_gene_snp_pairs_for_plotting.txt"), quote=F, sep="\t", col.names = F)
temp1 <- temp[which(temp$celltype=="KRT"),] %>% dplyr::select(gene, snp)
fwrite(temp1, paste0(dir.out,"/temp_KRT_gene_snp_pairs_for_plotting.txt"), quote=F, sep="\t", col.names = F)

# go to each celltype folder to plot each gene
# use /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/script/analysis2_fewIFNBeGene_investigation_plotGeneSnpPair.sh
# genes are plotted in /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/fewIFNBeGene_investigation/plots
# organize the individual gene:snp plots based on gene annotation
library(qpdf)

library(qpdf)

plot_dir  <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/fewIFNBeGene_investigation/plots"
merge_dir <- file.path(plot_dir, "merged")
dir.create(merge_dir, showWarnings = FALSE, recursive = TRUE)

temp <- left_join(
  tree_df_annotated,
  total_res[, c("snp", "gene", "celltype")]
)

# Keep track of folders we create
folders_seen <- c()

for (i in 1:nrow(temp)) {
  
  this_gene <- temp$gene[i]
  this_snp <- temp$snp[i]
  this_celltype <- temp$celltype[i]
  this_leaf_node <- temp$leaf_node[i]
  this_fail_reason <- temp$fail_reason[i]
  this_IFNB_eGene_status <- temp$is_IFNB_eGene[i]
  
  # Make correct folder name for this gene
  if (this_IFNB_eGene_status == FALSE) {
    this_folder_name <- paste0(
      "IFNG_induced_eGene_node",
      this_leaf_node,
      "_NOT_IFNB_eGene_",
      this_fail_reason
    )
  } else {
    this_folder_name <- paste0(
      "IFNG_induced_eGene_node",
      this_leaf_node,
      "_",
      this_fail_reason
    )
  }
  
  # Clean folder name
  this_folder_name <- sub(" \\(.*$", "", this_folder_name)
  this_folder_name <- sub("/.*$", "", this_folder_name)
  this_folder_name <- gsub(" ", "_", this_folder_name)
  
  folder_path <- file.path(plot_dir, this_folder_name)
  dir.create(folder_path, showWarnings = FALSE, recursive = TRUE)
  
  folders_seen <- unique(c(folders_seen, this_folder_name))
  
  # Find this gene's PDF
  pdf_path <- file.path(
    plot_dir,
    paste0(this_celltype, "_plot_", this_gene, "_", this_snp, ".pdf")
  )
  
  # Copy PDF into correct folder
  if (file.exists(pdf_path)) {
    file.copy(
      from = pdf_path,
      to = folder_path,
      overwrite = TRUE
    )
    
    cat("Copied:", basename(pdf_path), "->", this_folder_name, "\n")
  } else {
    cat("Missing PDF:", pdf_path, "\n")
  }
}

# Merge PDFs inside each folder
for (this_folder_name in folders_seen) {
  
  folder_path <- file.path(plot_dir, this_folder_name)
  
  pdfs_in_folder <- list.files(
    folder_path,
    pattern = "\\.pdf$",
    full.names = TRUE
  )
  
  if (length(pdfs_in_folder) > 0) {
    output_pdf <- file.path(
      merge_dir,
      paste0(this_folder_name, "_merged.pdf")
    )
    
    pdf_combine(
      input = pdfs_in_folder,
      output = output_pdf
    )
    
    cat("Merged", length(pdfs_in_folder), "files ->", output_pdf, "\n")
  } else {
    cat("No PDFs to merge in:", this_folder_name, "\n")
  }
}

# node 2 genes that are not IFNB eGenes ----
target_genes <- temp %>%
  dplyr::filter(leaf_node==2 & is_IFNB_eGene==F) %>%
  dplyr::select(c("gene","snp","celltype"))

# are those genes induced by IFNB? log2FC_IFNB or binary T/F
ifnb_lfc <- deseq_res_all %>%
  dplyr::filter(condition == "IFNB") %>%
  dplyr::select(gene, celltype, log2FC_IFNB = log2FoldChange, padj_IFNB = padj)

ifng_lfc <- deseq_res_all %>%
  dplyr::filter(condition == "IFNG") %>%
  dplyr::select(gene, celltype, log2FC_IFNG = log2FoldChange, padj_IFNG = padj)

target_genes %>%
  left_join(ifnb_lfc, by = c("gene", "celltype")) %>%
  left_join(ifng_lfc, by = c("gene", "celltype")) %>%
  group_by(celltype) %>%
  summarise(
    n               = n(),
    pct_induced_fc0.5     = mean(log2FC_IFNB > 0.5, na.rm = TRUE),
    pct_induced_fc1 = mean(log2FC_IFNB > 1,   na.rm = TRUE),
    median_log2FC   = median(log2FC_IFNB,      na.rm = TRUE)
  )

# Do we see a hint of eQTL? same direction beta?  Is beta_IFNB same sign as beta_IFNG? What's the beta_IFNB / beta_IFNG ratio? A concordance rate well above 50% would suggest real shared signal being masked by power.
target_genes2 <- left_join(target_genes, ifnb_lfc) %>%
  left_join(. , ifng_lfc) %>%
  left_join(. , total_res) %>%
  mutate(
    same_direction = sign(eQTL_beta_IFNG) == sign(eQTL_beta_IFNB),
    ifnb_induced = log2FC_IFNB > 1
  ) 

View(target_genes2 %>%
  group_by(celltype, ifnb_induced) %>%
  summarise(
    n              = n(),
    pct_same_dir   = mean(same_direction, na.rm = TRUE),
    median_eQTL_beta_ratio_IFNB_over_IFNG = median(abs(eQTL_beta_IFNB) / abs(eQTL_beta_IFNG), na.rm = TRUE),
    median_log2FC_ratio_IFNB_over_IFNG = median(log2FC_IFNB / log2FC_IFNG)
  )
)

# genes that aren't induced by IFNB but have same eQTL beta concordance -- what are they?
View(target_genes2 %>%
  dplyr::filter(ifnb_induced == FALSE & same_direction == TRUE)
)

# Controls — induced in similar magnitude but no eQTL in IFNG, similar donor number


##### ----- IFNB low eGene issue: decision tree modeling [>>> anchor separately on IFNG, IFNB, TNF eGenes <<<] ----
dir.out   <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/fewIFNBeGene_investigation/"
dir.in    <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap/"
dir.rna   <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/"
dir.model <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/"

celltypes  <- c("FRB", "KRT", "MEL")
cytokines  <- c("IFNG", "IFNB", "TNF")

# LOAD SHARED DATA ----
metagene_dict <- data.table::fread("/pi/manuel.garber-umw/human/skin/eQTLs/literature/metaidname.txt") %>% 
  dplyr::select(-id) %>% distinct()
metagene_dict_collapsed <- metagene_dict %>% group_by(meta) %>% 
  summarize(source=paste(name,collapse=','))
convert_meta_genes <- function(gene_vector) {
  # Identify meta genes in the input vector
  meta_genes <- gene_vector[grepl("meta", gene_vector)]
  non_meta_genes <- gene_vector[!grepl("meta", gene_vector)]
  
  # Convert meta genes to actual gene names using metagene_dict
  meta_to_gene <- metagene_dict$name[metagene_dict$meta %in% meta_genes]
  
  # Combine non-meta genes and converted meta genes, ensuring uniqueness
  converted_genes <- unique(c(non_meta_genes, meta_to_gene))
  
  return(converted_genes)
}

# big eGene annotation table
df <- read_tsv(paste0(dir.in, "/cisQTL_sigSNPs_ATACpeak_exprDEinduced_annotated.txt"))

# DEGs
load(paste0(dir.rna, "all_degs_abslog2FC1_padj0.05_post_outlier_exclusion.RData"))
deg_all <- bind_rows(
  deg.frb.ifnb %>% mutate(celltype="FRB", condition="IFNB"),
  deg.frb.ifng %>% mutate(celltype="FRB", condition="IFNG"),
  deg.frb.tnf  %>% mutate(celltype="FRB", condition="TNF"),
  deg.krt.ifnb %>% mutate(celltype="KRT", condition="IFNB"),
  deg.krt.ifng %>% mutate(celltype="KRT", condition="IFNG"),
  deg.krt.tnf  %>% mutate(celltype="KRT", condition="TNF"),
  deg.mel.ifnb %>% mutate(celltype="MEL", condition="IFNB"),
  deg.mel.ifng %>% mutate(celltype="MEL", condition="IFNG"),
  deg.mel.tnf  %>% mutate(celltype="MEL", condition="TNF")
)
rm(deg.frb.ifnb, deg.frb.ifng, deg.frb.tnf,
   deg.krt.ifnb, deg.krt.ifng, deg.krt.tnf,
   deg.mel.ifnb, deg.mel.ifng, deg.mel.tnf)

# DESeq2 results
load(paste0(dir.rna, "all_DESeq2_results_post_outlier_exclusion.RData"))
deseq_res_all <- bind_rows(
  as.data.frame(res.frb.ifnb) %>% rownames_to_column("gene") %>% mutate(celltype="FRB", condition="IFNB"),
  as.data.frame(res.frb.ifng) %>% rownames_to_column("gene") %>% mutate(celltype="FRB", condition="IFNG"),
  as.data.frame(res.frb.tnf)  %>% rownames_to_column("gene") %>% mutate(celltype="FRB", condition="TNF"),
  as.data.frame(res.krt.ifnb) %>% rownames_to_column("gene") %>% mutate(celltype="KRT", condition="IFNB"),
  as.data.frame(res.krt.ifng) %>% rownames_to_column("gene") %>% mutate(celltype="KRT", condition="IFNG"),
  as.data.frame(res.krt.tnf)  %>% rownames_to_column("gene") %>% mutate(celltype="KRT", condition="TNF"),
  as.data.frame(res.mel.ifnb) %>% rownames_to_column("gene") %>% mutate(celltype="MEL", condition="IFNB"),
  as.data.frame(res.mel.ifng) %>% rownames_to_column("gene") %>% mutate(celltype="MEL", condition="IFNG"),
  as.data.frame(res.mel.tnf)  %>% rownames_to_column("gene") %>% mutate(celltype="MEL", condition="TNF")
)
rm(res.frb.ifnb, res.frb.ifng, res.frb.tnf,
   res.krt.ifnb, res.krt.ifng, res.krt.tnf,
   res.mel.ifnb, res.mel.ifng, res.mel.tnf)

# CPM
cpm_matrix <- fread(paste0(dir.rna, "CPM.sampleFiltered.metaConverted.txt")) %>%
  dplyr::select(-name, -gene) %>%
  dplyr::filter(final_gene != "ARMCX5-GPRASP2") %>%
  column_to_rownames("final_gene")

meta_matrix <- fread(paste0(dir.rna, "metadata.sampleFiltered.txt"))

cpm_long <- cpm_matrix %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to="sample", values_to="cpm") %>%
  left_join(meta_matrix %>% select(sample, celltype, condition), by="sample") %>%
  group_by(gene, celltype, condition) %>%
  summarise(meanCPM = mean(cpm), .groups="drop") %>%
  pivot_wider(names_from=condition, values_from=meanCPM, names_prefix="meanCPM_")

# phyloP
# in bash, run this:
# workingdir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/fewIFNBeGene_investigation
# SNPbed=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap/QTL_sig_snps_all_celltype_conditions.bed
# prefix=sig_snps_all_celltype_conditions
# bash /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/05_get_phylop.sh ${workingdir} ${SNPbed} ${prefix}
phyloP_table <- fread(paste0(dir.out, "/phyloP_scores_sig_snps_all_celltype_conditions.txt"),
                      header=FALSE) %>%
  set_colnames(c("snp","size","bp_covered","phyloP_sum","phyloP_mean_snp","phyloP_mean_covered")) %>%
  dplyr::select(snp, phyloP = phyloP_mean_covered)

# eGene sets per cytokine x celltype
egene_sets <- map_dfr(celltypes, function(ct) {
  map_dfr(cytokines, function(cyt) {
    df %>%
      filter(celltype==ct, condition==cyt, QTLtype=="eQTL") %>%
      distinct(gene) %>%
      mutate(celltype=ct, anchor_condition=cyt)
  })
})

# induced gene sets per cytokine x celltype
induced_sets <- map_dfr(celltypes, function(ct) {
  map_dfr(cytokines, function(cyt) {
    deg_all %>%
      filter(celltype==ct, condition==cyt, padj<0.05, log2FoldChange>1.5) %>%
      distinct(gene) %>%
      mutate(celltype=ct, condition=cyt)
  })
})

# BUILD EXPANDED FEATURE MATRIX (all 3 anchors) ----
total_res_expanded <- data.frame()

for (anchor_cyt in cytokines) {
  message("\n========== Anchor: ", anchor_cyt, " ==========")
  
  for (ct in celltypes) {
    
    # anchor: induced eGenes for this cytokine x celltype
    anchor_induced <- induced_sets %>%
      filter(celltype==ct, condition==anchor_cyt) %>%
      pull(gene) %>%
      convert_meta_genes()
    
    anchor_egenes <- egene_sets %>%
      filter(celltype==ct, anchor_condition==anchor_cyt) %>%
      pull(gene) %>%
      convert_meta_genes()
    
    anchor_induced_egenes <- intersect(anchor_induced, anchor_egenes)
    
    # eGene sets for other cytokines (for is_*_eGene flags)
    ifnb_egenes <- egene_sets %>% filter(celltype==ct, anchor_condition=="IFNB") %>% pull(gene)
    ifng_egenes <- egene_sets %>% filter(celltype==ct, anchor_condition=="IFNG") %>% pull(gene)
    tnf_egenes  <- egene_sets %>% filter(celltype==ct, anchor_condition=="TNF")  %>% pull(gene)
    
    message(ct, ": ", length(anchor_induced_egenes), " ", anchor_cyt, " induced eGenes")
    
    if (length(anchor_induced_egenes) == 0) next
    
    # chunk file
    chunk_file <- ifelse(
      ct == "MEL",
      paste0(dir.model, ct, "/data/gene_chunk_dict_old.txt"),
      paste0(dir.model, ct, "/data/gene_chunk_dict.txt")
    )
    
    chunkF <- fread(chunk_file) %>%
      filter(gene %in% anchor_induced_egenes) %>%
      mutate(chunk = sprintf("%03d", chunk))
    
    res_this <- data.frame()
    
    for (this_gene in chunkF$gene) {
      
      # lead SNP = anchor cytokine's lead SNP for this gene
      this_snp <- df %>%
        filter(celltype==ct, condition==anchor_cyt, QTLtype=="eQTL", gene==this_gene) %>%
        arrange(p) %>%
        dplyr::slice(1) %>%
        pull(snp)
      
      if (length(this_snp) == 0) next
      
      chunk <- chunkF %>% 
        filter(gene==this_gene) %>% 
        dplyr::slice(1) %>% 
        pull(chunk)
      
      this_res <- fread(paste0(dir.model, ct, "/results/result_", chunk, ".tsv")) %>%
        filter(gene==this_gene, snp==this_snp) %>%
        dplyr::select(
          snp, gene,
          eQTL_beta_IFNG, eQTL_p_IFNG,
          eQTL_beta_IFNB, eQTL_p_IFNB,
          eQTL_beta_TNF,   eQTL_p_TNF,
          qc_n_per_condition_by_genotype
        ) %>%
        mutate(celltype=ct, anchor=anchor_cyt) %>%
        tidyr::separate(
          qc_n_per_condition_by_genotype,
          into = c(
            "g0_PBS","g0_IFNG","g0_IFNB","g0_TNF",
            "g1_PBS","g1_IFNG","g1_IFNB","g1_TNF",
            "g2_PBS","g2_IFNG","g2_IFNB","g2_TNF"
          ),
          sep="_", convert=TRUE
        ) %>%
        dplyr::select(-contains("PBS")) %>%
        left_join(
          df %>% dplyr::select(-tag, -condition, -QTLtype, -beta, -p),
          by=c("snp","gene","celltype")
        )
      
      res_this <- bind_rows(res_this, this_res)
    }
    
    # add eGene status flags
    res_this <- res_this %>%
      mutate(
        is_IFNB_eGene = gene %in% ifnb_egenes,
        is_IFNG_eGene = gene %in% ifng_egenes,
        is_TNF_eGene  = gene %in% tnf_egenes
      ) %>%
      distinct()
    
    total_res_expanded <- bind_rows(total_res_expanded, res_this)
  }
}


# JOIN CPM + phyloP ----

total_res_expanded <- total_res_expanded %>%
  left_join(cpm_long, by=c("gene","celltype")) %>%
  mutate(
    logCPM_PBS  = log2(meanCPM_PBS  + 1),
    logCPM_IFNB = log2(meanCPM_IFNB + 1),
    logCPM_IFNG = log2(meanCPM_IFNG + 1),
    logCPM_TNF  = log2(meanCPM_TNF  + 1)
  ) %>%
  left_join(phyloP_table, by="snp")

# HANDLE DUPLICATES ----
# a gene can appear in multiple anchors — keep all rows, anchor column tracks origin
# check overlap
total_res_expanded %>%
  dplyr::count(anchor, celltype) %>%
  print()

total_res_expanded %>%
  dplyr::count(gene, celltype) %>%
  filter(n > 1) %>%
  nrow() %>%
  message("Genes appearing in >1 anchor: ", .)

# WRITE OUTPUT ----
fwrite(total_res_expanded,
       file  = paste0(dir.out, "/all3anchors_inducedEGenes_eQTL_DESeq2_ATAC_phyloP_feature_matrix.txt"),
       sep   = "\t",
       quote = FALSE)

message("Done. Rows: ", nrow(total_res_expanded),
        " | Unique gene x celltype: ",
        n_distinct(paste(total_res_expanded$gene, total_res_expanded$celltype)))

message("is_IFNB_eGene distribution:")
print(table(total_res_expanded$is_IFNB_eGene, total_res_expanded$anchor))

# EXPANDED DECISION TREE (all 3 anchors, bio features only) ----


total_res_expanded <- fread(paste0(dir.out, "/all3anchors_inducedEGenes_eQTL_DESeq2_ATAC_phyloP_feature_matrix.txt"))

# prep feature matrix ----
tree_df_expanded <- total_res_expanded %>%
  dplyr::rename(log2FoldChange_FRB_IFNG = log2FoldChange_PBS_IFNG) %>%
  rowwise() %>%
  mutate(
    # gene induction
    induced_IFNB = get(paste0("is_induced_", celltype, "_IFNB")),
    induced_IFNG = get(paste0("is_induced_", celltype, "_IFNG")),
    induced_TNF  = get(paste0("is_induced_", celltype, "_TNF")),
    # ATAC openness — celltype matched
    open_IFNG        = get(paste0("open_", celltype, "_IFNG")),
    open_PBS         = get(paste0("open_", celltype, "_PBS")),
    peak_log2FC_IFNG = get(paste0("log2FoldChange_", celltype, "_IFNG")),
    peak_padj_IFNG   = get(paste0("padj_",           celltype, "_IFNG")),
    peak_dynamic     = get(paste0("peakDynamic_",     celltype, "_IFNG"))
  ) %>%
  ungroup() %>%
  left_join(
    deseq_res_all %>%
      filter(condition == "IFNB") %>%
      dplyr::select(gene, celltype, log2FC_IFNB = log2FoldChange),
    by = c("gene", "celltype")
  ) %>%
  left_join(
    deseq_res_all %>%
      filter(condition == "IFNG") %>%
      dplyr::select(gene, celltype, log2FC_IFNG = log2FoldChange),
    by = c("gene", "celltype")
  ) %>%
  left_join(
    deseq_res_all %>%
      filter(condition == "TNF") %>%
      dplyr::select(gene, celltype, log2FC_TNF = log2FoldChange),
    by = c("gene", "celltype")
  ) %>%
  mutate(
    # QTL eGene status
    is_IFNB_eGene = as.factor(is_IFNB_eGene),
    is_IFNG_eGene = as.factor(is_IFNG_eGene),   
    is_TNF_eGene  = as.factor(is_TNF_eGene),
    # gene induction
    induced_IFNB  = as.factor(as.character(induced_IFNB)),
    induced_IFNG  = as.factor(as.character(induced_IFNG)),
    induced_TNF   = as.factor(as.character(induced_TNF)),
    celltype      = as.factor(celltype),
    log2FC_IFNB   = replace_na(log2FC_IFNB, 0),
    log2FC_IFNG   = replace_na(log2FC_IFNG, 0),
    log2FC_TNF    = replace_na(log2FC_TNF,  0),
    # derived power features
    g2_diff_IFNG_IFNB = g2_IFNG - g2_IFNB,
    g2_diff_TNF_IFNB  = g2_TNF  - g2_IFNB,
    # ATAC: NA means no overlapping peak → treat as FALSE/0
    has_ATACpeak     = !is.na(ATACpeak_chr),
    open_IFNG        = ifelse(is.na(open_IFNG),        FALSE, open_IFNG),
    open_PBS         = ifelse(is.na(open_PBS),          FALSE, open_PBS),
    peak_log2FC_IFNG = ifelse(is.na(peak_log2FC_IFNG),  0,     peak_log2FC_IFNG),
    peak_padj_IFNG   = ifelse(is.na(peak_padj_IFNG),    1,     peak_padj_IFNG),
    peak_dynamic     = ifelse(is.na(peak_dynamic), "no_peak", peak_dynamic),
    # convert ATAC logical/character to factor
    has_ATACpeak  = as.factor(has_ATACpeak),
    open_IFNG     = as.factor(as.character(open_IFNG)),
    open_PBS      = as.factor(as.character(open_PBS)),
    peak_dynamic  = as.factor(peak_dynamic),
    # phyloP: impute median instead of dropping rows
    phyloP = ifelse(is.na(phyloP), median(phyloP, na.rm=TRUE), phyloP)
  ) %>%
  dplyr::select(
    is_IFNB_eGene,
    # genotype power
    g2_IFNG, g2_IFNB, g2_TNF,
    g2_diff_IFNG_IFNB, g2_diff_TNF_IFNB,
    # transcription
    induced_IFNB, induced_IFNG, induced_TNF,
    log2FC_IFNB, log2FC_IFNG, log2FC_TNF,
    # expression
    logCPM_PBS, logCPM_IFNB, logCPM_IFNG, logCPM_TNF,
    # conservation
    phyloP,
    # anchor eQTL effect
    eQTL_beta_IFNG, eQTL_beta_TNF,
    # eGene status for other cytokines
    is_IFNG_eGene, is_TNF_eGene,
    # ATAC
    has_ATACpeak,
    open_IFNG, open_PBS,
    peak_log2FC_IFNG, peak_padj_IFNG,
    peak_dynamic,
    # celltype
    celltype
  ) %>%
  drop_na()

# check balance after drop_na
message("Rows after drop_na: ", nrow(tree_df_expanded))
print(table(tree_df_expanded$is_IFNB_eGene))

# check how many NA were dropped (mostly from phyloP?)
message("Rows dropped: ", nrow(total_res_expanded) - nrow(tree_df_expanded))

# fit tree ----
set.seed(42)
tree_fit_expanded <- rpart(
  is_IFNB_eGene ~ .,
  data    = tree_df_expanded,
  method  = "class",
  control = rpart.control(
    maxdepth  = 4,
    minsplit  = 10,
    minbucket = 5,
    cp        = 0.01
  ),
  parms = list(prior = c(
    "FALSE" = sum(tree_df_expanded$is_IFNB_eGene == "FALSE") / nrow(tree_df_expanded),
    "TRUE"  = sum(tree_df_expanded$is_IFNB_eGene == "TRUE")  / nrow(tree_df_expanded)
  ))
)

# plot tree ----
pdf(paste0(dir.out, "/decision_tree_expanded_all3anchors_with_ATAC.pdf"), width=16, height=9)
rpart.plot(
  tree_fit_expanded,
  type        = 4,
  extra       = 104,
  box.palette = list("Purples", "Greens"),
  shadow.col  = "grey80",
  nn          = TRUE,
  main        = "Decision tree (all 3 anchors + ATAC): what predicts IFNB eGene status?"
)
dev.off()

# feature importance ----
tree_fit_expanded$variable.importance %>%
  as.data.frame() %>%
  rownames_to_column("feature") %>%
  set_colnames(c("feature", "importance")) %>%
  arrange(desc(importance)) %>%
  ggplot(aes(x = reorder(feature, importance), y = importance)) +
  geom_col(fill = "#4472C4") +
  coord_flip() +
  labs(title = "Feature importance (all 3 anchors + ATAC)",
       x = NULL, y = "Importance") +
  theme_classic(base_size = 12)

ggsave(paste0(dir.out, "/decision_tree_expanded_feature_importance_with_ATAC.pdf"),
       width = 7, height = 7)

# print tree structure ----
print(tree_fit_expanded)

# annotate genes based on expanded decision tree leaf nodes ----

tree_df_expanded$leaf_node <- as.integer(
  rownames(tree_fit_expanded$frame)[tree_fit_expanded$where]
)

node_label_map <- c(
  # FALSE leaves
  "8"  = "NOT_IFNB_eGene: low IFNB induction + low power + weak IFNG effect",
  "18" = "NOT_IFNB_eGene: low IFNB induction + low power + strong IFNG effect + low TNF expr",
  "10" = "NOT_IFNB_eGene: low IFNB induction + adequate power + weak IFNG effect",
  "12" = "NOT_IFNB_eGene: IFNB induced but TNF eGene competition",
  "26" = "NOT_IFNB_eGene: IFNB+TNF co-induced, not TNF eGene, low SNP conservation",
  "28" = "NOT_IFNB_eGene: IFNB induced + low TNF power + weak TNF effect",
  # TRUE leaves
  "19" = "IS_IFNB_eGene: rescued by strong IFNG effect + high TNF expression",
  "11" = "IS_IFNB_eGene: rescued by strong IFNG effect + adequate power",
  "27" = "IS_IFNB_eGene: IFNB+TNF co-induced, conserved SNP",
  "29" = "IS_IFNB_eGene: IFNB induced + strong negative TNF effect",
  "15" = "IS_IFNB_eGene: IFNB induced + adequate TNF power"
)

tree_df_expanded_annotated <- tree_df_expanded %>%
  mutate(
    gene      = total_res_expanded$gene,
    snp       = total_res_expanded$snp,
    anchor    = total_res_expanded$anchor,
    node_label = node_label_map[as.character(leaf_node)]
  )

# summary table ----
tree_df_expanded_annotated %>%
  dplyr::count(is_IFNB_eGene, node_label, celltype) %>%
  arrange(is_IFNB_eGene, node_label, celltype) %>%
  print(n = 50)

# bar chart ----
false_colors <- c(
  "NOT_IFNB_eGene: low IFNB induction + low power + weak IFNG effect"              = "#D0CFEE",
  "NOT_IFNB_eGene: low IFNB induction + low power + strong IFNG effect + low TNF expr" = "#A09FD8",
  "NOT_IFNB_eGene: low IFNB induction + adequate power + weak IFNG effect"         = "#7F77DD",
  "NOT_IFNB_eGene: IFNB induced but TNF eGene competition"                         = "#C0392B",
  "NOT_IFNB_eGene: IFNB+TNF co-induced, not TNF eGene, low SNP conservation"       = "#E8A090",
  "NOT_IFNB_eGene: IFNB induced + low TNF power + weak TNF effect"                 = "#F4A460"
)

true_colors <- c(
  "IS_IFNB_eGene: rescued by strong IFNG effect + high TNF expression"             = "#A8D8A8",
  "IS_IFNB_eGene: rescued by strong IFNG effect + adequate power"                  = "#5DBB63",
  "IS_IFNB_eGene: IFNB+TNF co-induced, conserved SNP"                             = "#1D9E75",
  "IS_IFNB_eGene: IFNB induced + strong negative TNF effect"                       = "#117A65",
  "IS_IFNB_eGene: IFNB induced + adequate TNF power"                               = "#0E6655"
)

all_colors <- c(false_colors, true_colors)

tree_df_expanded_annotated %>%
  dplyr::count(node_label, celltype) %>%
  mutate(
    outcome    = if_else(startsWith(node_label, "IS"), "IS IFNB eGene", "NOT IFNB eGene"),
    node_label = factor(node_label, levels = names(all_colors))
  ) %>%
  ggplot(aes(x = n, y = node_label, fill = node_label)) +
  geom_col() +
  geom_text(aes(label = n), hjust = -0.2, size = 3) +
  scale_fill_manual(values = all_colors) +
  facet_wrap(~celltype, nrow = 1) +
  labs(
    title = "Gene counts by decision tree node — expanded 3-anchor model",
    x = "n genes", y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position  = "none",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    axis.text.y      = element_text(size = 9)
  ) +
  xlim(0, max(
    tree_df_expanded_annotated %>% dplyr::count(node_label, celltype) %>% pull(n)
  ) * 1.15)

ggsave(paste0(dir.out, "/expanded_tree_node_annotation_barplot.pdf"), width = 14, height = 6)

# write annotated table ----
fwrite(
  tree_df_expanded_annotated %>%
    dplyr::select(gene, snp, celltype, anchor, is_IFNB_eGene, leaf_node, node_label),
  paste0(dir.out, "/expanded_tree_gene_node_annotations.txt"),
  sep = "\t", quote = FALSE
)

# write gene snp pairs for plotting ----
temp <- left_join(tree_df_expanded_annotated, total_res[,c("snp","gene","celltype")])

temp1 <- temp[which(temp$celltype=="MEL"),] %>% dplyr::select(gene, snp)
fwrite(temp1, paste0(dir.out,"/temp_MEL_gene_snp_pairs_for_plotting.txt"), quote=F, sep="\t", col.names = F)
temp1 <- temp[which(temp$celltype=="FRB"),] %>% dplyr::select(gene, snp)
fwrite(temp1, paste0(dir.out,"/temp_FRB_gene_snp_pairs_for_plotting.txt"), quote=F, sep="\t", col.names = F)
temp1 <- temp[which(temp$celltype=="KRT"),] %>% dplyr::select(gene, snp)
fwrite(temp1, paste0(dir.out,"/temp_KRT_gene_snp_pairs_for_plotting.txt"), quote=F, sep="\t", col.names = F)

# go to each celltype folder to plot each gene
# use /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/script/analysis2_fewIFNBeGene_investigation_plotGeneSnpPair.sh
# genes are plotted in /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/fewIFNBeGene_investigation/plots
# organize the individual gene:snp plots based on gene annotation
library(qpdf)

# 1. clean folder names
node_folder_map <- c(
  # FALSE leaves
  "8"  = "n8_NOT_IFNB_eGene_lowInduction_lowPower_weakIFNG",
  "18" = "n18_NOT_IFNB_eGene_lowInduction_lowPower_strongIFNG_lowTNFexpr",
  "10" = "n10_NOT_IFNB_eGene_lowInduction_adequatePower_weakIFNG",
  "12" = "n12_NOT_IFNB_eGene_IFNBinduced_TNFeGene_competition",
  "26" = "n26_NOT_IFNB_eGene_IFNBTNFcoinduced_notTNFeGene_lowPhyloP",
  "28" = "n28_NOT_IFNB_eGene_IFNBinduced_lowTNFpower_weakTNFeffect",
  # TRUE leaves
  "19" = "n19_IS_IFNB_eGene_rescued_strongIFNG_highTNFexpr",
  "11" = "n11_IS_IFNB_eGene_rescued_strongIFNG_adequatePower",
  "27" = "n27_IS_IFNB_eGene_IFNBTNFcoinduced_conservedSNP",
  "29" = "n29_IS_IFNB_eGene_IFNBinduced_strongNegTNFeffect",
  "15" = "n15_IS_IFNB_eGene_IFNBinduced_adequateTNFpower"
)

plot_dir  <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/fewIFNBeGene_investigation/plots"
merge_dir <- file.path(plot_dir, "merged")
dir.create(merge_dir, showWarnings = FALSE)

# 2. for each gene, find its PDF, copy to folder, then merge
temp %>%
  mutate(folder = node_folder_map[as.character(leaf_node)]) %>%
  group_by(folder) %>%
  group_walk(function(df, grp) {
    
    folder_path <- file.path(plot_dir, grp$folder)
    dir.create(folder_path, showWarnings = FALSE)
    
    # copy PDFs
    pdfs <- file.path(plot_dir, paste0(df$celltype, "_plot_", df$gene, "_", df$snp, ".pdf"))
    existing <- pdfs[file.exists(pdfs)]
    file.copy(existing, folder_path)
    
    # merge
    if (length(existing) > 0) {
      pdf_combine(input  = existing,
                  output = file.path(merge_dir, paste0(grp$folder, "_merged.pdf")))
      cat("Merged", length(existing), "files →", grp$folder, "\n")
    }
  })

# phyloP by induction specificity ----


phyloP_df <- total_res_expanded %>%
  filter(!is.na(phyloP)) %>%
  rowwise() %>%
  mutate(
    induced_IFNB = get(paste0("is_induced_", celltype, "_IFNB")),
    induced_IFNG = get(paste0("is_induced_", celltype, "_IFNG")),
    induced_TNF  = get(paste0("is_induced_", celltype, "_TNF"))
  ) %>%
  ungroup() %>%
  mutate(
    n_conditions_induced = as.integer(induced_IFNB) +
      as.integer(induced_IFNG) +
      as.integer(induced_TNF),
    induction_pattern = case_when(
      induced_IFNB & !induced_IFNG & !induced_TNF ~ "IFNB only",
      induced_IFNG & !induced_IFNB & !induced_TNF ~ "IFNG only",
      induced_TNF  & !induced_IFNB & !induced_IFNG ~ "TNF only",
      induced_IFNB &  induced_IFNG & !induced_TNF ~ "IFNB + IFNG",
      induced_IFNB &  induced_TNF  & !induced_IFNG ~ "IFNB + TNF",
      induced_IFNG &  induced_TNF  & !induced_IFNB ~ "IFNG + TNF",
      induced_IFNB &  induced_IFNG &  induced_TNF  ~ "All three",
      TRUE ~ "none"
    ),
    induction_pattern = factor(induction_pattern, levels = c(
      "IFNB only", "IFNG only", "TNF only",
      "IFNB + IFNG", "IFNB + TNF", "IFNG + TNF", "All three", "none"
    ))
  )

# summary table
phyloP_df %>%
  group_by(induction_pattern) %>%
  summarise(
    n             = n(),
    mean_phyloP   = round(mean(phyloP),   3),
    median_phyloP = round(median(phyloP), 3),
    pct_positive  = round(mean(phyloP > 0) * 100, 1),
    pct_above2    = round(mean(phyloP > 2) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(median_phyloP)) %>%
  print()

# wilcoxon tests: each pattern vs "none" as reference
patterns <- levels(phyloP_df$induction_pattern)
patterns <- patterns[patterns != "none"]

map_dfr(patterns, function(pat) {
  x <- phyloP_df %>% filter(induction_pattern == pat) %>% pull(phyloP)
  y <- phyloP_df %>% filter(induction_pattern == "none") %>% pull(phyloP)
  if (length(x) < 3) return(NULL)
  tibble(
    pattern = pat,
    n = length(x),
    median_phyloP = median(x),
    p_vs_none = wilcox.test(x, y)$p.value
  )
}) %>%
  arrange(p_vs_none) %>%
  print()

# violin plot ----
phyloP_df %>%
  filter(induction_pattern != "none") %>%
  ggplot(aes(x = induction_pattern, y = phyloP, fill = induction_pattern)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 2, linetype = "dotted", color = "grey70") +
  scale_fill_manual(values = c(
    "IFNB only"  = "#4472C4",
    "IFNG only"  = "#ED7D31",
    "TNF only"   = "#70AD47",
    "IFNB + IFNG"= "#9B59B6",
    "IFNB + TNF" = "#3498DB",
    "IFNG + TNF" = "#E67E22",
    "All three"  = "#C0392B"
  )) +
  facet_wrap(~celltype, nrow = 1) +
  labs(
    title    = "phyloP conservation score by induction specificity",
    subtitle = "Dashed = neutral (0); dotted = meaningful conservation (2)",
    x        = NULL, y = "phyloP20way score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 35, hjust = 1),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold")
  )

ggsave(paste0(dir.out, "/phyloP_by_induction_specificity.pdf"), width = 12, height = 5)

# also: single vs shared induction (collapsed) ----
phyloP_df %>%
  filter(induction_pattern != "none") %>%
  mutate(
    specificity = case_when(
      n_conditions_induced == 1 ~ "single condition",
      n_conditions_induced == 2 ~ "two conditions",
      n_conditions_induced == 3 ~ "all three conditions"
    )
  ) %>%
  ggplot(aes(x = specificity, y = phyloP, fill = specificity)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c(
    "single condition"      = "#7F77DD",
    "two conditions"        = "#1D9E75",
    "all three conditions"  = "#C0392B"
  )) +
  facet_wrap(~celltype, nrow = 1) +
  stat_compare_means(comparisons = list(
    c("single condition", "two conditions"),
    c("single condition", "all three conditions"),
    c("two conditions",   "all three conditions")
  ), method = "wilcox.test", label = "p.signif") +
  labs(
    title = "phyloP by number of conditions inducing the gene",
    x = NULL, y = "phyloP20way score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position  = "none",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold")
  )

ggsave(paste0(dir.out, "/phyloP_by_induction_breadth.pdf"), width = 10, height = 5)






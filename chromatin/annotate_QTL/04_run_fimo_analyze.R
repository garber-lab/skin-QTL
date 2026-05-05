library(dplyr)
library(readr)
library(stringr)
path.Rscript.fetch_expressed_genes="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_expressed_genes.R"
source(path.Rscript.fetch_expressed_genes)
path.Rscript.fetch_DE_genes="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_DE_genes_log2FC1_padj0.05.R"
source(path.Rscript.fetch_DE_genes)
path.Rscript.fetch_induced_genes="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_induced_genes_log2FC1_padj0.05.R"
source(path.Rscript.fetch_induced_genes)
trans_acting_genes_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/literature/trans_acting_genes/compiled_trans_acting_candidate_genes_and_category.txt"
trans_acting_genes_table <- fread(trans_acting_genes_file) %>%
  set_colnames(c("gene","gene_category"))
path.Rscript.fetch_mean_CPM_of_genes="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_mean_CPM_of_genes.R"
source(path.Rscript.fetch_mean_CPM_of_genes)

# -----------------------------
# user settings
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
this_snp <- args[2]
this_celltype <- args[3] 
this_condition <- args[4]

dir <- "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/rs838146"
this_snp <- "rs838146"
this_celltype <- "FRB"
this_condition <- "TNF"

fimo_file <- paste0(dir, "/fimo_output_", this_snp, ".txt")
out_file <- paste0(dir, "/fimo_output_analysis_",this_snp,"_",this_celltype,"_",this_condition,".txt")

snp_seq_ref <- fread(paste0(dir,"/SNP_slop100_REF.fa"))
snp_seq_alt <- fread(paste0(dir,"/SNP_slop100_ALT.fa"))

# -----------------------------
# read FIMO comparison result
# -----------------------------
fimo <- read_tsv(fimo_file, show_col_types = FALSE) %>%
  dplyr::mutate(TF = str_replace(motif_id, "\\.H12CORE.*$", ""))

# -----------------------------
# extract TF gene symbol from motif_id
# examples:
#   EGR1.H12CORE.0.PS.A  -> EGR1
#   SP1.H12CORE.0.P.B    -> SP1
# -----------------------------
lambert_tf <- fread("/pi/manuel.garber-umw/human/skin/eQTLs/literature/Lambert_2018_human_TF_table.csv", sep=",") %>%
  dplyr::select(c("Name","DBD")) %>%
  set_colnames(c("Name","TF_family"))
lambert_tf$TF_family <- gsub(";_", ";", gsub(" ", "_", lambert_tf$TF_family))

fimo2 <- fimo %>%
  left_join( . , lambert_tf, by=c("TF"="Name")) %>%
  dplyr::filter(
    # meaningful allelic difference (~5 fold change in p-value)
    priority_abs_delta_log10p >= 0.75
  )

# -----------------------------
# summarize to TF level
# one TF can have multiple motif models / rows
# -----------------------------
tf_summary <- fimo2 %>%
  group_by(TF, TF_family) %>%
  summarise(
    n_motifs = n(),
    n_gained_in_ALT = sum(motif_change_class == "gained_in_ALT", na.rm = TRUE),
    n_lost_in_ALT = sum(motif_change_class == "lost_in_ALT", na.rm = TRUE),
    n_stronger_in_ALT = sum(motif_change_class == "stronger_in_ALT", na.rm = TRUE),
    n_weaker_in_ALT = sum(motif_change_class == "weaker_in_ALT", na.rm = TRUE),
    best_priority_abs_delta_log10p = max(priority_abs_delta_log10p, na.rm = TRUE),
    best_abs_delta_score = max(abs(delta_score_ALT_minus_REF), na.rm = TRUE),
    strongest_change_class = motif_change_class[which.max(replace_na(priority_abs_delta_log10p, -Inf))],
    best_qvalue = min(pmin(qvalue_REF, qvalue_ALT), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    best_priority_abs_delta_log10p = ifelse(is.infinite(best_priority_abs_delta_log10p), NA, best_priority_abs_delta_log10p),
    best_abs_delta_score = ifelse(is.infinite(best_abs_delta_score), NA, best_abs_delta_score)
  )

# -----------------------------
gene_list <- unique(tf_summary$TF)

ct <- this_celltype
cond <- this_condition
expressed_genes <- unique(fetch_expressed_genes(ct, cond))
meanCPM <- fetch_mean_CPM_of_genes(ct, cond) %>%
  dplyr::filter(gene %in% gene_list)

if (cond == "PBS") {
  DE_genes <- character(0)
  induced_genes <- character(0)
} else {
  DE_genes <- unique(fetch_DE_genes(ct, cond))
  induced_genes <- unique(fetch_induced_genes(ct, cond))
}

gene_exprs_annotation <- data.frame(
  gene = gene_list,
  celltype = ct,
  condition = cond,
  is_expressed = gene_list %in% expressed_genes,
  is_DE = if (cond == "PBS") FALSE else gene_list %in% DE_genes,
  is_induced = if (cond == "PBS") FALSE else gene_list %in% induced_genes
) %>%
  left_join( . , meanCPM, by="gene")

# -----------------------------
# annotate with expression / DE / induction
# -----------------------------
tf_ranked <- tf_summary %>%
  left_join(gene_exprs_annotation, by = c("TF"="gene")) %>%
  mutate(
     interesting_by_motif = case_when(
      n_gained_in_ALT > 0 ~ TRUE,
      n_lost_in_ALT > 0 ~ TRUE,
      best_priority_abs_delta_log10p >= 0.5 ~ TRUE,
      best_abs_delta_score >= 3 ~ TRUE,
      TRUE ~ FALSE
    ),
    interesting_biologically = is_expressed,
    interesting_final = interesting_by_motif & interesting_biologically
  ) %>%
  arrange(
    desc(interesting_final),
    desc(n_gained_in_ALT + n_lost_in_ALT),
    desc(best_priority_abs_delta_log10p),
    desc(best_abs_delta_score)
  ) %>%
  mutate(
    priority_tier = case_when(
      interesting_final & best_priority_abs_delta_log10p >= 1.0 ~ "Tier1_high_confidence",
      interesting_final & best_priority_abs_delta_log10p >= 0.75 ~ "Tier2_moderate",
      interesting_final ~ "Tier3_expressed_hits",
      TRUE ~ "not_prioritized"
    )
  )
tf_ranked$priority_tier <- ordered(tf_ranked$priority_tier, levels=c("Tier1_high_confidence","Tier2_moderate","Tier3_expressed_hits","not_prioritized"))

# -----------------------------
# write a full ranked table
# -----------------------------
tf_ranked_filtered <- tf_ranked %>%
  dplyr::filter(interesting_final==TRUE) %>%
  arrange(priority_tier)


fwrite(tf_ranked_filtered, file=out_file, quote=F, sep="\t")
cat("Wrote fimo TF result to:", out_file, "\n")

# -------------------------------
# save RDS for downstream processing
# -------------------------------
out_prefix <- paste0(dir, "/fimo_summary_", this_snp, "_", this_celltype, "_", this_condition)
save(fimo, fimo2, tf_ranked, file = paste0(out_prefix, "objects.RData"))

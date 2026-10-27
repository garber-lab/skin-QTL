#!/usr/bin/env Rscript
Sys.setenv(TZ = "America/New_York")
suppressPackageStartupMessages({ library(data.table); library(stringr) })

#### --------- PER-CONTEXT CHROMATIN-DRIVER SCREEN ------------ ####

## ---- params ----
args      <- commandArgs(trailingOnly = TRUE)
this_gene <- args[1] 
#this_gene <- "AASS"
letter    <- substr(this_gene, 1, 1)
workingdir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats/", letter, "/", this_gene)

infile  <- file.path(workingdir, "big_table.txt")
out_summary <- file.path(workingdir, "mechanistic_variant_table_detail.txt")
out_candidate_list <- file.path(workingdir, "mechanistic_variant_list.txt")

## ---- read ----
dt <- fread(infile, na.strings = c("NA", ""))

## ---- helper function ----
add_allele_orientation <- function(big_table) {
  stopifnot(
    all(c("REF", "ALT", "cbp_allele1", "cbp_allele2") %in% colnames(big_table))
  )
  
  dt <- data.table::copy(big_table)
  
  dt[, ref_len := nchar(REF)]
  dt[, alt_len := nchar(ALT)]
  
  dt[, allele_orientation := data.table::fcase(
    # SNPs: direct allele comparison
    ref_len == 1 & alt_len == 1 & REF == cbp_allele1,
    "same",
    
    ref_len == 1 & alt_len == 1 & REF == cbp_allele2,
    "flipped",
    
    # Insertions in SNP genotype: REF shorter than ALT
    # Concordant if REF maps to "-" in chromBPNet, i.e. cbp_allele1 == "-"
    ref_len < alt_len & cbp_allele1 == "-",
    "same",
    
    ref_len < alt_len & cbp_allele2 == "-",
    "flipped",
    
    # Deletions in SNP genotype: REF longer than ALT
    # Concordant if ALT maps to "-" in chromBPNet, i.e. cbp_allele2 == "-"
    ref_len > alt_len & cbp_allele2 == "-",
    "same",
    
    ref_len > alt_len & cbp_allele1 == "-",
    "flipped",
    
    default = "mismatch_or_unhandled"
  )]
  
  dt[, c("ref_len", "alt_len") := NULL]
  
  dt
}

## ---- function to find candidate for a given context (celltype x condition x QTLtype) ----
# two parameters to fine-tune in step 3 and 4 of this function:
# 1. cbp_sig_threshold < 0.05 (chromBPNet metric pvalue threshold); 
# 2. acc_quantile_threshold <- 0.75 (chromBPNet active allele quantile being more than 75% of all peaks)
get_candidate_table_for_context <- function(this_celltype, this_condition, this_QTLtype) {
  
  #this_celltype="KRT"; this_condition="IFNB"; this_QTLtype="eQTL"
  this_dt <- dt %>%
    dplyr::filter(celltype == this_celltype, 
                  condition == this_condition, 
                  QTLtype == this_QTLtype)
  if (nrow(this_dt)==0) {return(NULL)}
  
  # 1. eGene? skip this context if not.
  if (unique(this_dt$is_egene)==FALSE) {return(NULL)}
  
  # 2. Filter sig QTLs.
  this_dt <- this_dt %>% dplyr::filter(this_dt$is_sig_snp==TRUE)
  if (nrow(this_dt)==0) {return(NULL)}
  
  # 3. Does the variant change chromatin?
  cbp_sig_threshold <- 0.05    # soft flag for "significant" chromatin disruption (NOT a hard filter)
  this_dt_filtered <- this_dt %>%
    dplyr::filter(this_dt$abs_logfc_x_active_allele_quantile.mean.pval < cbp_sig_threshold |
                    this_dt$jsd_x_active_allele_quantile.mean.pval < cbp_sig_threshold)
  if (nrow(this_dt_filtered)==0) {return(NULL)}
  
  # 4. Is this variant in a regulatory locus, existing or potential?
  #    OR logic: 
  #       peak open in THIS context             OR  
  #       peak gets opened up from PBS          OR
  #       active allele predicted accessible 
  #    The active-allele term is what captures peak CREATION (e.g. REF peak closed, ALT opens it).
  acc_quantile_threshold <- 0.75   # "high" = top quartile of genome-wide accessibility; tune
  this_dt_filtered_2 <- this_dt_filtered %>%
    dplyr::mutate(
      accessibility_evidence = 
        (ATACpeak_isopen %in% TRUE) |
        (!is.na(ATACpeak_padj) &  ATACpeak_padj < 0.05 & ATACpeak_log2FC > 0) |
        (!is.na(active_allele_quantile.mean) & active_allele_quantile.mean >= acc_quantile_threshold) 
    )
  
  # 5. Concordance? Is the accessibility-increasing allele = the expression-increasing allele?
  #    (concordant means the accessibility-increasing allele is the expression-increasing one)
  this_dt_filtered_3 <- add_allele_orientation(this_dt_filtered_2)
  this_dt_filtered_3[, cbp_sign := fcase(allele_orientation == "same",  1,
                                         allele_orientation == "flipped",    -1,
                                         default = NA_real_)]
  this_dt_filtered_3[, logfc_alt := logfc.mean * cbp_sign]            # now on REF->ALT axis, like beta
  this_dt_filtered_3[, cbp_eqtl_concordant := sign(logfc_alt) == sign(beta)]
  
  # 6. Rank the candidates by aggregating, in this priority order: 
  #     chromatin-effect strength (best of the two pvals) first, 
  #     accessibility (peak overlap as binary + active_allele_quantile as continuous) second, eQTL p third.
  #     best chromBPNet evidence by either axis (lower pval = stronger)
  this_dt_filtered_3[, cbp_best_pval := pmin(`abs_logfc_x_active_allele_quantile.mean.pval`,
                                             `jsd_x_active_allele_quantile.mean.pval`, na.rm = TRUE)]
  this_dt_filtered_3[, priority_rank := order(order(
    cbp_best_pval,                          # 1) strongest chromatin disruption first
    -as.integer(accessibility_evidence),    # 2) accessibility-supported first...
    -active_allele_quantile.mean,           #    ...then more accessible
    p,                                       # 3) strongest eQTL as the tiebreak
    na.last = TRUE))]
  setorder(this_dt_filtered_3, priority_rank)
  
  # Return ordered table
  candidate_table <- this_dt_filtered_3
  return(candidate_table)
}

## ---- loop through contexts ----
main_candidate_table <- data.frame()
for (ct in c("MEL", "KRT", "FRB")) {
  for (cond in c("PBS", "IFNG", "IFNB", "TNF")) {
    for (QTLtype in c("eQTL", "reQTL")) {
      if (cond=="PBS" & QTLtype=="reQTL"){next}
      
      #ct="FRB";cond="PBS";QTLtype="eQTL"
      this_candidate_table <- get_candidate_table_for_context(ct, cond, QTLtype)
      if (!is.null(this_candidate_table)) {
        main_candidate_table <- rbind(main_candidate_table, this_candidate_table)
      }

    }
  }
}

## ---- write output files ----
if (nrow(main_candidate_table) > 0) {
  fwrite(main_candidate_table, file = out_summary, sep = "\t", quote = FALSE)
  candidate_list <- unique(main_candidate_table$snp)
  writeLines(candidate_list, con = out_candidate_list)
} else {
  cat("", file = out_summary)
  cat("", file = out_candidate_list)
}
#message(paste0("wrote output to: ", out_summary))
#message(paste0("wrote output to: ", out_candidate_list))
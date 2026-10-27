# performs Benjamini Hochberg multiple testing correction across all genomic loci using eigenMT locally corrected p.values
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(magrittr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("
Usage:
  Rscript step5_BHcorrection.R <infile>

Example:
  Rscript adaptive_perm_gene_eQTL.R KRT_IFNB_eQTL.eigenMT.txt", "\n")
  quit(save = "no", status = 1)
}

infile         <- args[[1]]
basename <- unlist(strsplit(infile, ".eigenMT.txt"))
ct <- unlist(strsplit(basename, "_"))[1]
condition  <- unlist(strsplit(basename, "_"))[2]
QTLtype    <- unlist(strsplit(basename, "_"))[3]

dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct,"/eigenMT/results")

# read file 
f <- fread(paste0(dir,"/",infile))

# performs BH multiple testing correction
f[, p_gene_eigenMT := -expm1(Meff * log1p(-pmin))]
f$q_gene = p.adjust(f$p_gene_eigenMT, method="BH")
f = f %>% arrange(q_gene,pmin)
f_gene_fdr05 <- f[which(f$q_gene<0.05),] %>% arrange(q_gene,pmin)

# write results that passed fdr < 0.05
outname <- paste0(basename,"_gene_fdr05_table.txt")  
data.table::fwrite(f_gene_fdr05, 
                   file = paste0(dir,"/",outname), 
                   quote=F, sep="\t")

write(f_gene_fdr05$gene, file=paste0(dir,"/",basename,"_gene_fdr05_genelist.txt"))

# overwrite eigenMT table by adding the BH-corrected q-value back
data.table::fwrite(
  f,
  file = paste0(dir, "/", infile),
  quote = FALSE,
  sep = "\t"
)

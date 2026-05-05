suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
this_snp <- args[2]
output_file <- paste0(dir,"/SNP_context_insight.txt")

# dir <- "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/rs838146"
# this_snp <- "rs838146"
# output_file <- paste0(dir,"/SNP_context_insight.txt")

snp_seq_ref=fread(paste0(dir,"/SNP_slop100_REF.fa"))
snp_seq_alt=fread(paste0(dir,"/SNP_slop100_ALT.fa"))

# Extract SNP position from header
header_ref <- names(snp_seq_ref)  # or colnames, depending on your data.table structure
snp_pos_in_seq <- as.integer(str_extract(header_ref, "(?<=snp_pos=)\\d+"))

# Extract full sequences
ref_full <- strsplit(as.character(snp_seq_ref[1, 1]), "")[[1]]
alt_full <- strsplit(as.character(snp_seq_alt[1, 1]), "")[[1]]

# Get allele identities at SNP position
ref_allele <- ref_full[snp_pos_in_seq]
alt_allele <- alt_full[snp_pos_in_seq]

# ─── SEQUENCE INSIGHTS DATA.FRAME ────────────────────────────────────────────
local_window <- 10
local_ref <- unlist(ref_full)[(snp_pos_in_seq - local_window):(snp_pos_in_seq + local_window)]
local_alt <- unlist(alt_full)[(snp_pos_in_seq - local_window):(snp_pos_in_seq + local_window)]

# ── 1. Mutation type ──────────────────────────────────────────────────────────
purines     <- c("A","G")
pyrimidines <- c("C","T")
is_transition <- (ref_allele %in% purines     & alt_allele %in% purines) |
  (ref_allele %in% pyrimidines & alt_allele %in% pyrimidines)
mut_type <- ifelse(is_transition, "transition", "transversion")

# ── 2. GC content ─────────────────────────────────────────────────────────────
gc_pct <- function(b) round(mean(unlist(b) %in% c("G","C")) * 100, 1)
gc_ref_val <- gc_pct(local_ref)
gc_alt_val <- gc_pct(local_alt)
gc_delta   <- gc_alt_val - gc_ref_val
gc_direction <- case_when(
  gc_delta >= 3  ~ "increase",
  gc_delta <= -3 ~ "decrease",
  TRUE           ~ "no_change"
)

# ── 3. Homopolymer runs ───────────────────────────────────────────────────────
max_run_info <- function(bases) {
  bases <- unlist(bases)
  r   <- rle(bases)
  idx <- which.max(r$lengths)
  list(len=r$lengths[idx], base=r$values[idx])
}
run_ref <- max_run_info(local_ref)
run_alt <- max_run_info(local_alt)
run_effect <- case_when(
  run_alt$len > run_ref$len & run_alt$len >= 3 ~ "created",
  run_ref$len > run_alt$len & run_ref$len >= 3 ~ "disrupted",
  run_ref$len >= 3 & run_alt$len >= 3 &
    run_ref$len == run_alt$len               ~ "maintained",
  TRUE                                         ~ "none"
)
run_base      <- ifelse(run_effect == "created",  run_alt$base, run_ref$base)
run_length    <- ifelse(run_effect == "created",  run_alt$len,  run_ref$len)
run_detail    <- case_when(
  run_effect == "none"       ~ NA_character_,
  TRUE ~ paste0(run_effect, " run of ", run_length, " ", run_base, "'s")
)

# ── 4. CpG dinucleotides ──────────────────────────────────────────────────────
count_cpg <- function(bases) {
  seq_str <- paste(unlist(bases), collapse="")
  lengths(regmatches(seq_str, gregexpr("CG", seq_str)))
}
cpg_ref_n <- count_cpg(local_ref)
cpg_alt_n <- count_cpg(local_alt)
cpg_delta <- cpg_alt_n - cpg_ref_n
cpg_effect <- case_when(
  cpg_delta > 0  ~ "created",
  cpg_delta < 0  ~ "disrupted",
  cpg_ref_n == 0 ~ "none_in_either",
  TRUE           ~ "unchanged"
)
cpg_detail <- case_when(
  cpg_effect == "created"       ~ paste0("created ", abs(cpg_delta), " CpG site",
                                         ifelse(abs(cpg_delta)>1,"s","")),
  cpg_effect == "disrupted"     ~ paste0("disrupted ", abs(cpg_delta), " CpG site",
                                         ifelse(abs(cpg_delta)>1,"s","")),
  cpg_effect == "none_in_either"~ "no CpG in either allele",
  cpg_effect == "unchanged"     ~ paste0(cpg_ref_n, " CpG site",
                                         ifelse(cpg_ref_n>1,"s",""), " unchanged")
)

# ── 5. Palindrome context ─────────────────────────────────────────────────────
rev_comp <- function(bases) {
  bases <- unlist(bases)
  comp  <- c(A="T", T="A", G="C", C="G")
  rev(comp[bases])
}
is_palindrome <- function(bases) {
  bases <- unlist(bases)
  if (length(bases) %% 2 != 0) return(FALSE)
  mid <- length(bases) / 2
  all(bases[1:mid] == rev_comp(bases[(mid+1):length(bases)]))
}
pal_window    <- 6
pal_ref_bases <- unlist(ref_full)[(snp_pos_in_seq - pal_window):(snp_pos_in_seq + pal_window - 1)]
pal_alt_bases <- unlist(alt_full)[(snp_pos_in_seq - pal_window):(snp_pos_in_seq + pal_window - 1)]
pal_ref <- tryCatch(is_palindrome(pal_ref_bases), error=function(e) FALSE)
pal_alt <- tryCatch(is_palindrome(pal_alt_bases), error=function(e) FALSE)
palindrome_effect <- case_when(
  !pal_ref &  pal_alt ~ "created",
  pal_ref & !pal_alt ~ "disrupted",
  pal_ref &  pal_alt ~ "maintained",
  TRUE                ~ "none_in_either"
)
palindrome_detail <- case_when(
  palindrome_effect == "created"       ~ "creates palindromic sequence (favors homodimer TF binding)",
  palindrome_effect == "disrupted"     ~ "disrupts palindromic sequence",
  palindrome_effect == "maintained"    ~ "palindromic sequence maintained in both alleles",
  palindrome_effect == "none_in_either"~ "no palindrome in either allele"
)

# ── 6. Assemble data.frame ────────────────────────────────────────────────────
snp_insight_df <- tibble(
  snp                    = this_snp,
  ref_allele             = ref_allele,
  alt_allele             = alt_allele,
  mut_type               = mut_type,               # transition / transversion
  gc_ref_pct             = gc_ref_val,             # local GC% in REF
  gc_alt_pct             = gc_alt_val,             # local GC% in ALT
  gc_delta_pct           = gc_delta,               # ALT - REF
  gc_direction           = gc_direction,           # increase / decrease / no_change
  homopolymer_effect     = run_effect,             # created / disrupted / maintained / none
  homopolymer_detail     = run_detail,             # e.g. "created run of 5 G's"
  cpg_ref_count          = cpg_ref_n,             # # CpG in local REF window
  cpg_alt_count          = cpg_alt_n,             # # CpG in local ALT window
  cpg_delta              = cpg_delta,             # ALT - REF
  cpg_effect             = cpg_effect,            # created / disrupted / unchanged / none_in_either
  cpg_detail             = cpg_detail,
  palindrome_ref         = pal_ref,               # logical
  palindrome_alt         = pal_alt,               # logical
  palindrome_effect      = palindrome_effect,     # created / disrupted / maintained / none_in_either
  palindrome_detail      = palindrome_detail
)

fwrite(snp_insight_df, output_file, quote=F, sep="\t")

# # ── 8. Auto-generate subtitle from data.frame ─────────────────────────────────
# subtitle_parts <- c(
#   paste0(ref_allele, "→", alt_allele, " ", mut_type),
#   if (gc_direction != "no_change")
#     paste0(gc_direction, "s local GC content (", gc_ref_val, "% → ", gc_alt_val, "%)"),
#   if (!is.na(run_detail)) run_detail,
#   cpg_detail,
#   palindrome_detail
# )
# auto_subtitle <- paste(Filter(Negate(is.null), subtitle_parts), collapse=" · ")


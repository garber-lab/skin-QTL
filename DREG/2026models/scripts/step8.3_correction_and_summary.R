#!/usr/bin/env Rscript
# =============================================================================
# step8.3_correction_and_summary.R   (genome-wide run)
#
# STEP 3. BH-correct the conditional tests from step8.2, finalize the per-direction
# class, derive a gene-level story per comparison, and summarize.
#
# CORRECTION: one conditional test per gene per direction (query lead conditioned
# on the reference context's signal set) -> BH ACROSS GENES, per comparison, per
# direction (eigenMT not applicable: no per-gene multi-SNP scan here). BH is
# applied to the conditional p (p_cond) within the "test_needed" set.
#
# Per-direction final class (query lead Q vs reference context R; r2 = LD to nearest R signal):
#   same_signal_primary     r2>=0.8 to R's PRIMARY signal
#   same_signal_secondary   r2>=0.8 to a SECONDARY signal of R           [rank shift]
#   independent_strong      r2<0.5,  sig in R, survives conditioning (q_cond<FDR)  [clearly distinct variant]
#   independent_weak        0.5<=r2<0.8, sig in R, survives conditioning          [distinct, partial LD]
#   context_specific        r2<0.8,  NOT even marginally sig in R         [R-absent element]
#   tagged_complex          r2<0.8,  sig in R, killed by conditioning     [tagged by R's signals combined]
#   query_no_signal         query context has no discovered signal for this gene
#
# Gene story (priority): CONTEXT_SPECIFIC (a novel_in_ref in either direction) >
#   DISTINCT_CAUSAL (both distinct_independent) > RANK_SHIFT > SHARED_PRIMARY >
#   SINGLE_CONTEXT_eGene > AMBIGUOUS.
#
# Input : results/conditional_analysis/pairwise/matched_*_vs_*.txt (step8.2)
# Output: results/conditional_analysis/pairwise/corrected/<same>          (finalized per gene)
#         results/conditional_analysis/SUMMARY_per_comparison.txt         (category counts)
#         results/conditional_analysis/INTERESTING_findings.txt           (context-specific + distinct + rank-shift)
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
stamp <- function(...) message(sprintf("[%s] %s", format(Sys.time(),"%Y-%m-%d %H:%M:%S"), paste0(...)))

args <- commandArgs(trailingOnly=TRUE)
ROOT <- if (length(args)>=1) args[1] else "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models"
FDR     <- if (length(args)>=2) as.numeric(args[2]) else 0.05
MARG    <- if (length(args)>=3) as.numeric(args[3]) else 0.05
LD_WEAK <- if (length(args)>=4) as.numeric(args[4]) else 0.5   # r2 split for independent weak vs strong
PROX_BP <- if (length(args)>=5) as.numeric(args[5]) else 2000  # a1<->b1 proximity flag (soft; not a filter)
DIR  <- file.path(ROOT,"results","conditional_analysis")
PW   <- file.path(DIR,"pairwise"); COR <- file.path(PW,"corrected"); dir.create(COR, showWarnings=FALSE, recursive=TRUE)
GENO_FILE <- file.path(ROOT,"results","genotype_of_eGene_tested_variants.txt")
files <- list.files(PW, pattern="^matched_.*_vs_.*\\.txt$", full.names=TRUE)
if (!length(files)) stop("no matched_*_vs_*.txt in ", PW)

# ---- lead SNP positions (for the a1<->b1 physical-separation flag) ----------
leads <- unique(unlist(lapply(files, function(f){ x<-fread(f, select=c("a1","b1")); c(x$a1,x$b1) })))
leads <- leads[!is.na(leads) & leads!=""]
idf <- tempfile(); writeLines(leads, idf)
gp  <- fread(cmd=sprintf("awk -F'\t' 'NR==FNR{a[$1];next} FNR==1||($3 in a)' %s %s", shQuote(idf), shQuote(GENO_FILE)),
             select=c("CHROM","POS","ID")); unlink(idf)
POSMAP <- setNames(gp$POS, gp$ID); CHRMAP <- setNames(gp$CHROM, gp$ID)
sep_bp <- function(a,b) ifelse(is.na(a)|is.na(b)|!(a %in% names(POSMAP))|!(b %in% names(POSMAP))|CHRMAP[a]!=CHRMAP[b],
                               NA_integer_, abs(POSMAP[a]-POSMAP[b]))

# r2 tiers (query lead vs nearest reference signal):
#   >=0.8 same signal (primary or secondary); test_needed cases have r2<0.8:
#   present & survives conditioning -> independent_strong (r2<LD_WEAK) / independent_weak (LD_WEAK<=r2<0.8)
#   not even marginally present     -> context_specific (ref-absent element)
#   present but killed by conditioning -> tagged_complex
finalize_dir <- function(class_pre, r2, pmarg, qcond) {
  fifelse(class_pre=="no_query_signal", "query_no_signal",
   fifelse(class_pre=="same_primary", "same_signal_primary",
    fifelse(class_pre=="matches_ref_secondary", "same_signal_secondary",
     fifelse(is.na(pmarg), "untestable",
      fifelse(pmarg > MARG, "absent_in_ref",
       fifelse(!is.na(qcond) & qcond < FDR,
               fifelse(is.na(r2) | r2 < LD_WEAK, "independent_strong", "independent_weak"),
               "tagged_complex"))))))
}
# CONTEXT_SPECIFIC gene story is renamed by which axis is being compared:
#   across_celltypes -> CELLTYPE_SPECIFIC ;  across_conditions -> CONDITION_SPECIFIC
ctx_story <- function(mode) fifelse(mode=="across_celltypes","CELLTYPE_SPECIFIC","CONDITION_SPECIFIC")

summ <- list(); interesting <- list()
for (f in files) {
  d <- fread(f); if (!nrow(d)) next
  # BH within test_needed set, each direction
  d[, b1vsA_qcond := NA_real_][, a1vsB_qcond := NA_real_]
  iB <- which(d$b1vsA_class=="test_needed" & is.finite(d$b1vsA_pcond)); if (length(iB)) d$b1vsA_qcond[iB] <- p.adjust(d$b1vsA_pcond[iB],"BH")
  iA <- which(d$a1vsB_class=="test_needed" & is.finite(d$a1vsB_pcond)); if (length(iA)) d$a1vsB_qcond[iA] <- p.adjust(d$a1vsB_pcond[iA],"BH")
  d[, b1vsA_final := finalize_dir(b1vsA_class, b1vsA_r2, b1vsA_pmarg, b1vsA_qcond)]
  d[, a1vsB_final := finalize_dir(a1vsB_class, a1vsB_r2, a1vsB_pmarg, a1vsB_qcond)]

  # gene-level story (priority); CONTEXT_SPECIFIC relabeled by comparison axis
  d[, gene_story := fifelse(b1vsA_final=="absent_in_ref" | a1vsB_final=="absent_in_ref", ctx_story(mode),
      fifelse(grepl("^independent", b1vsA_final) & grepl("^independent", a1vsB_final), "DISTINCT_CAUSAL",
       fifelse(b1vsA_final=="same_signal_secondary" | a1vsB_final=="same_signal_secondary", "RANK_SHIFT",
        fifelse(b1vsA_final=="same_signal_primary" | a1vsB_final=="same_signal_primary", "SHARED_PRIMARY",
         fifelse(b1vsA_final=="query_no_signal" | a1vsB_final=="query_no_signal", "SINGLE_CONTEXT_eGene",
          "AMBIGUOUS")))))]
  # which context carries the extra/specific signal (absent in the other)
  d[, novel_in := fifelse(b1vsA_final=="absent_in_ref" & a1vsB_final=="absent_in_ref", "BOTH",
      fifelse(b1vsA_final=="absent_in_ref", paste0(ctxB,"(b1 absent in A)"),
       fifelse(a1vsB_final=="absent_in_ref", paste0(ctxA,"(a1 absent in B)"), NA_character_)))]
  # physical separation of the two leads + SOFT proximity flag (leads <=PROX_BP apart are
  # likely the same regulatory element -> weak evidence for distinct/specific; NOT filtered here)
  d[, a1b1_sep_bp := sep_bp(a1, b1)]
  d[, a1b1_proximal_flag := !is.na(a1b1_sep_bp) & a1b1_sep_bp <= PROX_BP]
  fwrite(d, file.path(COR, basename(f)), sep="\t")

  cmp <- sub("^matched_","",sub("\\.txt$","",basename(f)))
  st <- d[, .N, gene_story]; setnames(st,"N","n")
  wide <- dcast(st, .~gene_story, value.var="n"); wide[, `.`:=NULL]
  summ[[f]] <- cbind(data.table(comparison=cmp, mode=d$mode[1], held=d$held[1], n_genes=nrow(d)), wide)

  hit <- d[gene_story %in% c("CELLTYPE_SPECIFIC","CONDITION_SPECIFIC","DISTINCT_CAUSAL","RANK_SHIFT"),
    .(comparison=cmp, mode, held, ctxA, ctxB, gene, gene_story, novel_in,
      a1, b1, same_lead, a1b1_sep_bp, a1b1_proximal_flag,
      b1vsA_final, b1vsA_rank, b1vsA_r2, b1vsA_pmarg, b1vsA_qcond,
      a1vsB_final, a1vsB_rank, a1vsB_r2, a1vsB_pmarg, a1vsB_qcond,
      nA_signals, nB_signals)]
  if (nrow(hit)) interesting[[f]] <- hit
}
S <- rbindlist(summ, use.names=TRUE, fill=TRUE); for (c in names(S)) if (is.numeric(S[[c]])) S[is.na(get(c)), (c):=0]
setorder(S, mode, held)
fwrite(S, file.path(DIR,"SUMMARY_per_comparison.txt"), sep="\t")
I <- rbindlist(interesting, use.names=TRUE, fill=TRUE)
if (nrow(I)) {
  I[, story_ord := match(gene_story, c("CELLTYPE_SPECIFIC","CONDITION_SPECIFIC","DISTINCT_CAUSAL","RANK_SHIFT"))]
  setorder(I, story_ord, comparison, gene); I[, story_ord := NULL]
}
fwrite(I, file.path(DIR,"INTERESTING_findings.txt"), sep="\t")

# ---- GENE-LEVEL collapse: one row per gene, its story labels + comparisons ---
if (nrow(I)) {
  GL <- I[, .(
    n_findings   = .N,
    stories      = paste(sort(unique(gene_story)), collapse=", "),
    comparisons  = paste(sprintf("%s[%s]", comparison, gene_story), collapse=", "),
    any_distinct = any(gene_story=="DISTINCT_CAUSAL"),
    any_proximal_flag = any(a1b1_proximal_flag %in% TRUE)   # a lead pair <=PROX_BP apart in some comparison
  ), by=gene][order(-any_distinct, -n_findings, gene)]
  fwrite(GL, file.path(DIR,"GENE_LEVEL_summary.txt"), sep="\t")
  stamp("wrote GENE_LEVEL_summary.txt (", nrow(GL), " genes)")
}

stamp("=== SUMMARY per comparison ==="); print(S)
stamp("=== interesting findings: ", nrow(I), " gene-comparisons ===")
if (nrow(I)) print(I[, .N, by=.(mode, gene_story)])
stamp("context-specific: ", if(nrow(I)) I[gene_story %in% c('CELLTYPE_SPECIFIC','CONDITION_SPECIFIC'), .N] else 0,
      " | distinct_causal: ", if(nrow(I)) I[gene_story=='DISTINCT_CAUSAL', .N] else 0,
      " | rank_shift: ", if(nrow(I)) I[gene_story=='RANK_SHIFT', .N] else 0)
stamp("DONE step8.3.")

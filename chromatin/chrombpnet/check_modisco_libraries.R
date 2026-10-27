suppressPackageStartupMessages({library(dplyr); library(stringr); library(purrr); library(tibble)})
CB <- "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet"
clean_tf <- function(x) toupper(sub("\\.(mouse|human)$", "", sub("_.*$", "", x), ignore.case = TRUE))
parse_motifs_html <- function(fn) {
  x <- paste(readLines(fn, warn = FALSE), collapse = "\n")
  rows  <- str_split(x, "<tr")[[1]][-1]
  cells <- map(rows, ~ str_trim(gsub("<[^>]+>", "", str_match_all(.x, "(?s)<t[hd][^>]*>(.*?)</t[hd]>")[[1]][, 2])))
  hdr <- cells[[1]]; body <- keep(cells[-1], ~ length(.x) == length(hdr))
  df <- as_tibble(do.call(rbind, body), .name_repair = "minimal"); names(df) <- hdr
  transmute(df, pattern, motif_tf = clean_tf(match0),
            qval0 = suppressWarnings(as.numeric(qval0)),
            num_seqlets = suppressWarnings(as.integer(num_seqlets)))
}
for (m in c("KRT_PBS","KRT_IFN","FRB_PBS","FRB_IFN","MEL_PBS","MEL_IFN")) {
  fn <- file.path(CB, sprintf("ATAC_%s_merged_alldonors", m), "biasFactor0.5",
                  "fold0/chrombpnet_model/evaluation/modisco_profile/motifs.html")
  if (!file.exists(fn)) { cat(sprintf("%-8s MISSING %s\n", m, fn)); next }
  d <- parse_motifs_html(fn)
  irf <- filter(d, str_detect(motif_tf, "^(IRF|STAT|PRDM1)"))
  cat(sprintf("%-8s patterns=%-4d  IRF/STAT patterns=%-3d  %s\n", m, nrow(d), nrow(irf),
              paste(sprintf("%s(n=%s)", irf$motif_tf, irf$num_seqlets), collapse = ", ")))
}

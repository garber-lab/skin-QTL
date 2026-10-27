#!/usr/bin/env bash
# Per cell type: one master peak list from the PBS + IFN MACS2 peaks.
# Summits within 250 bp are one element; keep the strongest summit (highest -log10 q) with its
# original narrowPeak row, so every peak keeps a real summit for the 2114 bp model window.
# Output: ATAC_<CT>_master_summits.narrowPeak = 10 narrowPeak columns + col 11 condition (PBS / IFN / both)

set -eo pipefail
module load bedtools

for CT in KRT FRB MEL; do
  TMP=$(mktemp)
  # 1. one summit per row from PBS + IFN; cols: 1-3 summit, 4-13 original narrowPeak row, 14 condition
  for C in PBS IFN; do
    awk -F'\t' -v c=$C 'BEGIN{OFS="\t"} $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {s = $2 + $10; print $1, s, s + 1, $0, c}' \
      ATAC_${CT}_${C}_merged_alldonors_peaks.narrowPeak
  done | sort -k1,1 -k2,2n > "${TMP}"

  # 2-4. group summits within 250 bp, keep the strongest (col 12 = -log10 q), label PBS / IFN / both
  bedtools cluster -d 250 -i "${TMP}" | awk -F'\t' 'BEGIN{OFS="\t"}
    { k = $15
      if (!(k in q) || $12 > q[k]) { q[k] = $12; row = $4; for (i = 5; i <= 13; i++) row = row OFS $i; best[k] = row }
      seen[k, $14] = 1 }
    END { for (k in best) {
            c = (seen[k, "PBS"] && seen[k, "IFN"]) ? "both" : (seen[k, "PBS"] ? "PBS" : "IFN")
            print best[k], c } }' | sort -k1,1 -k2,2n > ATAC_${CT}_master_summits.narrowPeak
  rm -f "${TMP}"

  echo "== ${CT}: $(wc -l < ATAC_${CT}_master_summits.narrowPeak) master peaks"
  cut -f11 ATAC_${CT}_master_summits.narrowPeak | sort | uniq -c
done

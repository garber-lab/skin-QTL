#!/bin/bash
#BSUB -J sort_index_hm
#BSUB -q short
#BSUB -W 4:00
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=6000]"
# ^ per-SLOT on this cluster: 4 slots x 6000 = 24 GB total. The previous 16000
#   asked for 64 GB to run an 8 GB sort, which only delays scheduling.
#BSUB -o /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/log/sort_index_hm_%J.out
#BSUB -e /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/log/sort_index_hm_%J.err
# ---------------------------------------------------------------------------
# Drop every unprefixed column that duplicates an hm_ column, coordinate-sort on
# hm_pos, bgzip, tabix.
#
# Result: 24 columns -> 14. Values are left exactly as deposited; only whole
# columns are removed. Nothing is recomputed, lifted over, or re-harmonised --
# EBI already did that, and the hm_ columns ARE its output.
#
# Submit with:   bsub < scripts/sort_index_hm.sh
# ---------------------------------------------------------------------------
set -euo pipefail

dir=${COLOC_DIR:-/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc}
trait=${TRAIT:-systemic_lupus_erythematosus_bentham2015}
IN=${IN:-26502338-GCST003156-EFO_0002690.h.tsv.gz}
OUT=${OUT:-${trait}.hg38.tsv.gz}

# Columns to delete, BY NAME: every unprefixed column that has an hm_ counterpart.
#
#   chromosome, base_pair_location          build37 -> use hm_chrom, hm_pos
#   other_allele, effect_allele             pre-swap  -> use hm_other_allele, hm_effect_allele
#   beta, odds_ratio, ci_lower, ci_upper    pre-flip  -> use hm_beta, hm_odds_ratio, hm_ci_*
#   effect_allele_frequency                 pre-swap  -> use hm_effect_allele_frequency
#   variant_id                              duplicate of hm_rsid
#
# What survives is every hm_* column plus exactly two unprefixed ones,
# standard_error and p_value. Those two have no hm_ counterpart because they
# do not need one: neither an allele swap nor a strand flip changes a standard
# error or a p-value. Only beta's sign changes. So after this, there is no
# column left in the file that can be read in the wrong orientation.
DROP="chromosome base_pair_location variant_id other_allele effect_allele beta odds_ratio ci_lower ci_upper effect_allele_frequency"

# bgzip/tabix are not on PATH on this cluster; they come from the htslib module.
module load htslib/1.16

cd "${dir}/${trait}"
mkdir -p "${dir}/log" .tmp
trap 'rm -rf .tmp' EXIT

# --- 1. drop the build37 columns -------------------------------------------
# Rows with hm_pos = NA are also removed: they have no GRCh38 coordinate, so
# they cannot be sorted or indexed. The count is reported below.
echo "[1/3] dropping columns: ${DROP}"
zcat "$IN" | awk -F'\t' -v OFS='\t' -v drop="$DROP" '
# CRLF tolerance. With Windows line endings awk (RS="\n") leaves a bare \r on the
# LAST field of every record, so column 24 reads as "effect_allele_frequency\r"
# and fails to match by name while all 23 others match. Strip it before anything
# else looks at the record.
{ sub(/\r$/, "") }
NR==1 {
  n = split(drop, d, " "); for (i = 1; i <= n; i++) kill[d[i]] = 1
  for (i = 1; i <= NF; i++) {             # normalise names: BOM, stray CR, stray whitespace
    h = $i; gsub(/\r/, "", h); sub(/^\xef\xbb\xbf/, "", h); gsub(/^[ \t]+|[ \t]+$/, "", h)
    name[i] = h
  }
  for (i = 1; i <= NF; i++) {
    if (name[i] in kill) { found[name[i]] = 1; continue }
    keep[++k] = i
    if (name[i] == "hm_pos") pcol = k     # position of hm_pos AFTER the drop
  }
  for (i = 1; i <= n; i++) if (!(d[i] in found)) { print "no such column: " d[i] > "/dev/stderr"; bad = 1 }
  if (!pcol) { print "no hm_pos column - is this really a harmonised file?" > "/dev/stderr"; bad = 1 }
  if (bad) {                              # dump what was actually parsed, so the next
    printf "header as parsed (%d fields):\n", NF > "/dev/stderr"   # mismatch diagnoses itself
    for (i = 1; i <= NF; i++) printf "  [%2d] |%s| len=%d\n", i, name[i], length(name[i]) > "/dev/stderr"
    exit 2
  }
  PCOL = keep[pcol]                       # hm_pos index in the INPUT record
  line = name[keep[1]]; for (i = 2; i <= k; i++) line = line OFS name[keep[i]]
  print line
  next
}
NF == 0 { next }
{
  if ($PCOL == "NA" || $PCOL == "") { na++; next }
  line = $(keep[1]); for (i = 2; i <= k; i++) line = line OFS $(keep[i])
  print line
  kept++
}
END {
  printf "      kept %d rows, dropped %d with hm_pos=NA\n", kept+0, na+0 > "/dev/stderr"
  if (kept + 0 == 0) { print "no rows survived - refusing to write an empty file" > "/dev/stderr"; exit 1 }
}' > .tmp/body.tsv

# --- 2. sort ----------------------------------------------------------------
# Resolve hm_chrom / hm_pos positions from the NEW header, so the sort keys are
# never stale if DROP changes.
CHR=$(head -1 .tmp/body.tsv | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="hm_chrom") print i}')
POS=$(head -1 .tmp/body.tsv | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="hm_pos")   print i}')
[ -n "$CHR" ] && [ -n "$POS" ] || { echo "cannot locate hm_chrom/hm_pos" >&2; exit 1; }
echo "[2/3] sorting on hm_chrom (col ${CHR}), hm_pos (col ${POS})"

{
  head -1 .tmp/body.tsv
  tail -n +2 .tmp/body.tsv \
    | LC_ALL=C sort -t $'\t' -k"${CHR},${CHR}V" -k"${POS},${POS}n" \
           -S 8G --parallel=4 -T .tmp
} | bgzip -@ 4 -c > "$OUT"

# --- 3. index ---------------------------------------------------------------
echo "[3/3] indexing"
tabix -f -s "$CHR" -b "$POS" -e "$POS" -S 1 "$OUT"

echo "done: $(pwd)/${OUT}"
zcat "$OUT" | head -1 | tr '\t' '\n' | nl
echo "rows: $(( $(zcat "$OUT" | wc -l) - 1 ))"
echo "query test (chr6 MHC): $(tabix "$OUT" 6:32000000-32100000 | wc -l) rows"

# ---------------------------------------------------------------------------
# WHY THE UNHARMONISED EFFECT COLUMNS GO TOO, not just the coordinates.
#
# rs4607995 in this deposit (hm_code 11) is the worked example:
#     other_allele=T   effect_allele=C   beta=-0.030459    <- as submitted
#     hm_other_allele=C hm_effect_allele=T hm_beta=+0.030459  <- harmonised
# Same magnitude, swapped alleles, opposite sign. Reading `beta` instead of
# `hm_beta` inverts the direction of effect for every variant EBI had to flip,
# and nothing errors -- so it is the same class of silent failure as reading
# base_pair_location instead of hm_pos, just in the effect rather than the
# coordinate.
#
# The 14 surviving columns:
#    1 hm_variant_id   2 hm_rsid        3 hm_chrom      4 hm_pos
#    5 hm_other_allele 6 hm_effect_allele  7 hm_beta    8 hm_odds_ratio
#    9 hm_ci_lower    10 hm_ci_upper   11 hm_effect_allele_frequency
#   12 hm_code        13 p_value       14 standard_error
#
# hm_effect_allele_frequency is present but NA for every variant in this
# study -- Bentham 2015 did not deposit allele frequencies. That is a property
# of the deposit, not of this script. coloc type="cc" with N and s does not
# need it; LD-based finemapping does, and would need a 1000G EUR / TOPMed
# annotation step first.
#
# To keep a dropped column, delete its name from DROP above. The tabix flags
# are derived from the post-drop header, so they follow automatically.
# ---------------------------------------------------------------------------


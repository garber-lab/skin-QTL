#!/bin/bash
#BSUB -J drop_N_ssc
#BSUB -q short
#BSUB -W 1:00
#BSUB -n 1
#BSUB -R "rusage[mem=4000]"
#BSUB -o /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/log/drop_N_ssc_%J.out
#BSUB -e /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/log/drop_N_ssc_%J.err
# ---------------------------------------------------------------------------
# Drop the N column from the Lopez-Isac 2019 SSc meta-analysis file.
#
# WHY. N is not a sample size. This file is PLINK --meta-analysis output
# (CHR BP SNP A1 A2 N P P(R) OR OR(R) Q I), where N is the number of the 14
# contributing GWAS cohorts in which the SNP was present. Measured over 250k
# variants it runs 2..14, and the floor of 2 is the paper's own rule --
# "we considered all the SNPs that were shared by at least two data sets"
# (Lopez-Isac 2019 Nat Commun, PMC6823490, Methods). Feeding it to coloc as a
# per-variant N would understate the study by three orders of magnitude.
#
# The real values for coloc, from the README:
#     N = 26,679   (9,095 cases + 17,584 controls)
#     s = 9095 / 26679 = 0.340867
#
# Input is space-PADDED, not tab-delimited. Output is proper TSV.
#
# Submit:  bsub < scripts/drop_N_column_ssc.sh
# or run directly -- it is one awk pass over 450 MB, well under a minute:
#          bash scripts/drop_N_column_ssc.sh
# ---------------------------------------------------------------------------
set -euo pipefail

dir=${COLOC_DIR:-/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc}
trait=${TRAIT:-systemic_sclerosis_lopez2019}
IN=${IN:-Lopez-Isac_prePMID_META_GWAS_SSc.meta.txt}
OUT=${OUT:-Lopez-Isac_META_GWAS_SSc.noN.tsv}

# Column(s) to delete, by name. To keep N but stop it reading as a sample size,
# set DROP="" and RENAME="N:n_cohorts" instead -- it is a usable QC filter
# (a variant seen in 2 of 14 cohorts is far weaker evidence than one seen in 14).
# ${DROP-N} not ${DROP:-N}: the colon form treats an empty string as unset, so
# DROP="" would still drop N and the RENAME route above could never be taken.
DROP=${DROP-N}
RENAME=${RENAME:-}

cd "${dir}/${trait}"
[ -s "$IN" ] || { echo "missing or empty input: $(pwd)/$IN" >&2; exit 1; }

awk -v OFS='\t' -v drop="$DROP" -v rename="$RENAME" '
NR==1 {
  n = split(drop, d, " ");   for (i = 1; i <= n; i++) kill[d[i]] = 1
  m = split(rename, r, " "); for (i = 1; i <= m; i++) { split(r[i], kv, ":"); newname[kv[1]] = kv[2] }
  for (i = 1; i <= NF; i++) {
    if ($i in kill) { found[$i] = 1; continue }
    keep[++k] = i
    hdr[k] = ($i in newname) ? newname[$i] : $i
  }
  for (i = 1; i <= n; i++) if (!(d[i] in found)) { print "no such column: " d[i] > "/dev/stderr"; bad = 1 }
  if (bad) {
    printf "header as parsed (%d fields):\n", NF > "/dev/stderr"
    for (i = 1; i <= NF; i++) printf "  [%2d] |%s|\n", i, $i > "/dev/stderr"
    exit 2
  }
  line = hdr[1]; for (i = 2; i <= k; i++) line = line OFS hdr[i]
  print line
  next
}
NF == 0 { next }
{
  line = $(keep[1]); for (i = 2; i <= k; i++) line = line OFS $(keep[i])
  print line
  rows++
}
END {
  # awk runs END even after `exit 2`, and an exit here would override that code
  # with 1 ("no rows") -- masking the real cause. Re-raise it first.
  if (bad) exit 2
  printf "kept %d data rows, %d columns\n", rows+0, k > "/dev/stderr"
  if (rows + 0 == 0) { print "no rows - refusing to write empty output" > "/dev/stderr"; exit 1 }
}' "$IN" > "$OUT"

echo "done: $(pwd)/${OUT}"
head -3 "$OUT"

# ---------------------------------------------------------------------------
# STILL OUTSTANDING before this file can go into coloc -- none of it done here:
#
#   A2  is '?' for 100% of variants (checked: 250k sampled, every N from 2 to
#       14, no exceptions). Not recoverable from this file. See notes below.
#   OR  needs beta = log(OR); there is no SE column, so se = |beta| / |z| with
#       z from P. Guard the OR==1 / P==1 case, where that is 0/0.
#       Use the pair (P, OR) or the pair (P(R), OR(R)) -- never one of each.
#   hg19 -> liftOver to hg38. rs3094315 sits at chr1:752566 here; hg38 is
#       chr1:817186.
# ---------------------------------------------------------------------------

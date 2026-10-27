#!/bin/bash
module load conda; source activate fastQTL
set -euo pipefail
shopt -s nullglob
d=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/motifscan/allhits_profile/MEL_IFN
mkdir -p ${d}

mh=(/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5/fold0/chrombpnet_model/auxiliary/interpret_subsample/*modisco_results_profile_scores.h5)
if [[ ${#mh[@]} -ne 1 ]]; then echo "ERROR: expected 1 MoDISco h5 for MEL_IFN, found ${#mh[@]}" >&2; exit 1; fi
h5s=(/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/motifscan/contribs/MEL_IFN/MEL_IFN_fold*.profile_scores.h5)
if [[ ${#h5s[@]} -eq 0 ]]; then
  echo "ERROR: no profile contributions for MEL_IFN (run step7_motifScan_contribs.sh first)" >&2; exit 1
fi
echo "MEL_IFN: ${#h5s[@]} fold(s), $(wc -l < /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5/data/downloads/peaks_no_blacklist.narrowPeak) peaks"; date

# peak regions with genomic coordinates, contributions averaged across folds (cached)
if [[ ! -s ${d}/regions.npz ]]; then
   /home/shuo.shan-umw/.conda/envs/finemo/bin/finemo extract-regions-chrombpnet-h5     -c "${h5s[@]}" -p /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_IFN_merged_alldonors/biasFactor0.5/data/downloads/peaks_no_blacklist.narrowPeak -w 1000 -o ${d}/regions.npz
fi
date

# every motif at once (no -I): same cost as calling a single motif
 /home/shuo.shan-umw/.conda/envs/finemo/bin/finemo call-hits -r ${d}/regions.npz -m ${mh[0]} -o ${d}   -N /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/motifscan/names_profile/MEL_IFN/names.txt -c 0.0001 -S 10000
echo "${mh[0]}" > ${d}/modisco_source.txt
date

# BED of instances + per-motif counts
awk -F'\t' 'NR>1 && $1!="NA" {print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$13}' ${d}/hits.tsv   | sort -k1,1 -k2,2n > ${d}/hits.bed
awk -F'\t' 'NR>1 {c[$6]++} END {for (k in c) print k"\t"c[k]}' ${d}/hits.tsv   | sort -k2,2nr > ${d}/motif_counts.tsv
echo "MEL_IFN: $(wc -l < ${d}/hits.bed) instances, $(wc -l < ${d}/motif_counts.tsv) motifs"
head -15 ${d}/motif_counts.tsv

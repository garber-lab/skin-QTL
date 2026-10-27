#!/bin/bash
# FIMO REF vs ALT motif scan around one variant (adapted from annotate_QTL/04_run_fimo.sh)
#
# Usage:  bash step6_fimo_snp.sh <variant_dir> [flank=50]
#   <variant_dir> = posterchild/<rsID>, must contain <rsID>.variant.tsv
#                   (chr, pos 1-based, REF, ALT, id; written by step6_run_snp_contrib.sh)
#
# Output: <variant_dir>/motifs/fimo/
#   <rsID>.fimo_REF_ALT.tsv   motifs overlapping the variant, REF vs ALT merged + classified
#                             (same columns as fimo_output_<snp>.txt from annotate_QTL, plus TF and
#                              pos_rel coordinates relative to the variant, 0 = variant)
#   <rsID>.fimo_all_hits.tsv  every hit with p < 1e-4 in the window, both alleles, pos_rel coords
#                             (context: lets hits be matched to Fi-NeMo spans later)
#   SNP_window_REF.fa / SNP_window_ALT.fa, variant_info.tsv
set -euo pipefail

vdir=$(readlink -f "$1")
flank=${2:-50}
rsid=$(basename "${vdir}")

############################ CONFIG ############################
CB=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet
DL=${CB}/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/data/downloads
GENOME_FA=${DL}/hg38.fa                  # same genome the ChromBPNet models were trained on
CHROMSIZES=${DL}/hg38.chrom.sizes
MEME_DB=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/motif_databases/HUMAN/HOCOMOCOv12_core_HUMAN_mono_meme_format.meme
PY=/home/shuo.shan-umw/miniconda3/envs/fastQTL/bin/python3   # has pandas/numpy
ALL_HITS_P=1e-4                          # p-value cutoff for the context table
# FIMO report from the annotate_QTL scripts (set RUN_REPORT=0 to skip)
RUN_REPORT=${RUN_REPORT:-1}
OLD_SCRIPTS=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL
R_SIF=/share/pkg/containers/rstudio_example/r_4.5.2.sif
R_LIBS=${HOME}/R/x86_64-pc-linux-gnu-library/4.2
REPORT_CONTEXTS="KRT:IFNG KRT:PBS FRB:IFNG FRB:PBS MEL:IFNG MEL:PBS"   # celltype:condition
################################################################
export PATH=/home/shuo.shan-umw/meme/bin:/home/shuo.shan-umw/meme/libexec/meme-5.5.5:$PATH
module load bedtools

varF=${vdir}/${rsid}.variant.tsv
[[ -s "${varF}" ]] || { echo "ERROR: missing ${varF}" >&2; exit 1; }
command -v fimo >/dev/null || { echo "ERROR: fimo not on PATH" >&2; exit 1; }

out=${vdir}/motifs/fimo
mkdir -p "${out}"
cd "${out}"
echo "[FIMO] ${rsid}  flank=+/-${flank} bp  db=$(basename ${MEME_DB})"; date

### 1. variant BED (0-based) and +/- flank window, then genome sequence
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$5,$3,$4}' "${varF}" > SNP.bed      # chr start end id ref alt
bedtools slop -b "${flank}" -i SNP.bed -g "${CHROMSIZES}" > SNP_window.bed
bedtools getfasta -fi "${GENOME_FA}" -bed SNP_window.bed -fo SNP_window.fa

### 2. REF / ALT sequences, with the same allele checks as 04_run_fimo.sh
"${PY}" - <<'PY'
import sys
chrom, start, end, snp_id, ref, alt = open("SNP.bed").readline().split()[:6]
start, end, ref, alt = int(start), int(end), ref.upper(), alt.upper()
if "," in alt:
    sys.exit(f"Multi-allelic variant not supported: ALT={alt}")
wchrom, wstart, wend = open("SNP_window.bed").readline().split()[:3]
wstart = int(wstart)
seq = "".join(l.strip() for l in open("SNP_window.fa") if not l.startswith(">")).upper()

idx0 = start - wstart            # 0-based variant index in the window
pos1 = idx0 + 1                  # 1-based (FIMO coordinates)
ref_span, alt_span = seq[idx0:idx0 + len(ref)], seq[idx0:idx0 + len(alt)]
if ref_span == ref:
    anchor_len = len(ref)
elif alt_span == alt:
    anchor_len = len(alt)
    print(f"[WARN] genome matches ALT, not REF, at {chrom}:{end}; keeping REF/ALT labels", file=sys.stderr)
else:
    sys.exit(f"Neither allele matches genome at {chrom}:{end}: REF={ref} ({ref_span}) ALT={alt} ({alt_span})")

left, right = seq[:idx0], seq[idx0 + anchor_len:]
seq_ref, seq_alt = left + ref + right, left + alt + right
vtype = ("SNV" if len(ref) == 1 and len(alt) == 1 else
         "deletion" if len(ref) > len(alt) else "insertion" if len(ref) < len(alt) else "MNV")

for allele, s in [("REF", seq_ref), ("ALT", seq_alt)]:
    # both keys: 04_annotate_SNP_context.R greps "snp_pos=", 04_run_fimo.sh wrote "variant_pos="
    header = f">{snp_id}|{allele}|snp_pos={pos1}|variant_pos={pos1}|variant_type={vtype}"
    for fn in (f"SNP_window_{allele}.fa", f"SNP_slop100_{allele}.fa"):   # second name: annotate_QTL scripts
        with open(fn, "w") as o:
            o.write(f"{header}\n{s}\n")
with open("snp_position_in_seq.txt", "w") as o:                          # read by 04_compile_fimo_output_plots.R
    o.write(str(pos1) + "\n")
with open("variant_interval_REF.txt", "w") as o: o.write(f"{pos1}\t{pos1 + len(ref) - 1}\n")
with open("variant_interval_ALT.txt", "w") as o: o.write(f"{pos1}\t{pos1 + len(alt) - 1}\n")
with open("variant_info.tsv", "w") as o:   # columns of annotate_QTL/04_run_fimo.sh, plus window info
    o.write("variant_id\tchrom\tstart\tend\tref\talt\tvariant_type\tvariant_pos1\tref_start1\tref_end1\t"
            "alt_start1\talt_end1\tref_seq_len\talt_seq_len\tpos\twindow_start0\twindow_len\n")
    o.write(f"{snp_id}\t{chrom}\t{start}\t{end}\t{ref}\t{alt}\t{vtype}\t{pos1}\t"
            f"{pos1}\t{pos1 + len(ref) - 1}\t{pos1}\t{pos1 + len(alt) - 1}\t"
            f"{len(seq_ref)}\t{len(seq_alt)}\t{end}\t{wstart}\t{len(seq)}\n")
print(f"[FIMO] {snp_id} {chrom}:{end} {ref}>{alt} ({vtype}), variant at window pos {pos1}/{len(seq)}")
PY

### 3. FIMO on REF and ALT, no threshold (filter afterwards, as before)
for allele in REF ALT; do
  fimo --oc "results_${allele}" --thresh 1 --verbosity 1 "${MEME_DB}" "SNP_window_${allele}.fa"
  awk 'BEGIN{FS=OFS="\t"} !/^#/ && NR>1 && $1!=""' "results_${allele}/fimo.tsv" > "fimo_all_${allele}.tsv"
done

### 4. merge REF/ALT, classify, add TF name + pos_rel coordinates
"${PY}" - "${rsid}" "${ALL_HITS_P}" <<'PY'
import sys
import numpy as np
import pandas as pd

rsid, all_p = sys.argv[1], float(sys.argv[2])
cols = ["motif_id", "motif_alt_id", "sequence_name", "start", "stop", "strand",
        "score", "pvalue", "qvalue", "matched_sequence"]
info = pd.read_csv("variant_info.tsv", sep="\t").iloc[0]
vpos1 = int(info.variant_pos1)

def read(allele):
    try:
        d = pd.read_csv(f"fimo_all_{allele}.tsv", sep="\t", header=None, names=cols)
    except pd.errors.EmptyDataError:
        d = pd.DataFrame(columns=cols)
    lo, hi = map(int, open(f"variant_interval_{allele}.txt").read().split())
    d["overlaps_variant"] = (d["start"] <= hi) & (d["stop"] >= lo)
    return d

ref, alt = read("REF"), read("ALT")

# ---- context table: all reasonably strong hits in the window, both alleles ----
ctx = pd.concat([ref.assign(allele="REF"), alt.assign(allele="ALT")])
ctx = ctx[ctx["pvalue"] < all_p].copy()
ctx["TF"] = ctx["motif_id"].str.replace(r"\.H12CORE.*$", "", regex=True)
ctx["start_rel"], ctx["stop_rel"] = ctx["start"] - vpos1, ctx["stop"] - vpos1
ctx.sort_values(["allele", "pvalue"]).to_csv(f"{rsid}.fimo_all_hits.tsv", sep="\t", index=False)

# ---- variant-overlapping hits: best per motif per allele, merged (as in 04_run_fimo.sh) ----
def best_per_motif(d, sfx):
    keep = ["start", "stop", "strand", "score", "pvalue", "qvalue", "matched_sequence"]
    d = d[d["overlaps_variant"]]
    if d.empty:
        return pd.DataFrame(columns=["motif_id"] + [f"{c}_{sfx}" for c in keep])
    d = d.sort_values(["motif_id", "pvalue"]).groupby("motif_id", as_index=False).first()
    return d[["motif_id"] + keep].rename(columns={c: f"{c}_{sfx}" for c in keep})

m = best_per_motif(ref, "REF").merge(best_per_motif(alt, "ALT"), on="motif_id", how="outer")
m["present_REF"] = m["score_REF"].notna()
m["present_ALT"] = m["score_ALT"].notna()

def classify(r):
    if r.present_REF and not r.present_ALT: return "lost_in_ALT"
    if r.present_ALT and not r.present_REF: return "gained_in_ALT"
    if r.present_REF and r.present_ALT:
        if r.score_ALT > r.score_REF: return "stronger_in_ALT"
        if r.score_ALT < r.score_REF: return "weaker_in_ALT"
        return "unchanged"
    return "absent_both"

m["neg_log10p_REF"] = -np.log10(m["pvalue_REF"].astype(float).clip(lower=1e-300))
m["neg_log10p_ALT"] = -np.log10(m["pvalue_ALT"].astype(float).clip(lower=1e-300))
m["delta_score_ALT_minus_REF"] = m["score_ALT"] - m["score_REF"]
m["delta_log10p_ALT_minus_REF"] = m["neg_log10p_ALT"] - m["neg_log10p_REF"]
m["priority_abs_delta_log10p"] = m["delta_log10p_ALT_minus_REF"].abs()
m["motif_change_class"] = m.apply(classify, axis=1)
m["TF"] = m["motif_id"].str.replace(r"\.H12CORE.*$", "", regex=True)
for a in ["REF", "ALT"]:
    m[f"start_rel_{a}"] = m[f"start_{a}"] - vpos1
    m[f"stop_rel_{a}"] = m[f"stop_{a}"] - vpos1

m = m.sort_values(["motif_change_class", "priority_abs_delta_log10p"],
                  ascending=[True, False], na_position="last")
m.insert(0, "variant_id", rsid)
m.to_csv(f"{rsid}.fimo_REF_ALT.tsv", sep="\t", index=False)

strong = m[m["priority_abs_delta_log10p"] >= 0.75]
print(f"[FIMO] {len(m)} motifs overlap the variant; {len(strong)} with |delta -log10p| >= 0.75")
if len(strong):
    print(strong[["TF", "motif_change_class", "pvalue_REF", "pvalue_ALT",
                  "delta_log10p_ALT_minus_REF", "start_rel_ALT", "stop_rel_ALT"]]
          .head(15).to_string(index=False))
PY

### 5. FIMO report PDFs from the annotate_QTL scripts, one set per celltype x condition
if [[ "${RUN_REPORT}" == "1" ]]; then
  cp "${rsid}.fimo_REF_ALT.tsv" "fimo_output_${rsid}.txt"    # name 04_run_fimo_analyze.R opens
  export R_LIBS_USER=${R_LIBS}
  singularity exec "${R_SIF}" Rscript "${OLD_SCRIPTS}/04_annotate_SNP_context.R" "${out}" "${rsid}" \
    || echo "[FIMO] WARNING: 04_annotate_SNP_context.R failed" >&2
  for ctx in ${REPORT_CONTEXTS}; do
    ct=${ctx%%:*}; cond=${ctx##*:}
    { singularity exec "${R_SIF}" Rscript "${OLD_SCRIPTS}/04_run_fimo_analyze.R" "${out}" "${rsid}" "${ct}" "${cond}" &&
      singularity exec "${R_SIF}" Rscript "${OLD_SCRIPTS}/04_compile_fimo_output_plots.R" "${out}" "${rsid}" "${ct}" "${cond}"; } \
      || echo "[FIMO] WARNING: report failed for ${ct} ${cond}" >&2
  done
fi

### 5b. copy the report PDFs next to the ChromBPNet figures, one name per context
if [[ "${RUN_REPORT}" == "1" ]]; then
  shopt -s nullglob
  n=0
  for f in fimo_summary_${rsid}_*.pdf; do
    ctx=${f#fimo_summary_${rsid}_}; ctx=${ctx%.pdf}          # e.g. KRT_IFNG
    cp -f "${f}" "${vdir}/${rsid}.fimo_report_${ctx}.pdf" && n=$((n+1))
  done
  echo "[FIMO] copied ${n} report PDF(s) to ${vdir}/${rsid}.fimo_report_<celltype>_<condition>.pdf"
fi

### 6. tidy up (keep final tables, fasta and variant_info for inspection)
rm -rf results_REF results_ALT fimo_all_REF.tsv fimo_all_ALT.tsv SNP_window.fa
echo "[FIMO] wrote ${out}/${rsid}.fimo_REF_ALT.tsv and ${rsid}.fimo_all_hits.tsv"; date

#!/bin/bash
set -euo pipefail
module load bedtools

workingdir=$1
snp_id=$2
snp_bed="${3:-}"

genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
genome_fa=/share/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.fa
memeF=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/motif_databases/HUMAN/HOCOMOCOv12_core_HUMAN_mono_meme_format.meme
export PATH=/home/shuo.shan-umw/meme/bin:/home/shuo.shan-umw/meme/libexec/meme-5.5.5:$PATH

cd "${workingdir}"
echo "fimo starting for snp" > "${workingdir}/fimo_status.txt"

########## create SNP bed if needed
if [[ -z "${snp_bed}" || ! -f "${snp_bed}" ]]; then
    echo "Warning: SNP bed not provided or not found, generating now"
    bash /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/01_compile_snp_bed.sh "${workingdir}" "${snp_id}"
    snp_bed="${workingdir}/QTL_${snp_id}.bed"
fi

### make 100bp flanking window around SNP: ±50 bp
bedtools slop -b 50 -i "${snp_bed}" -g "${genome}" > SNP_slop100.bed

### get fasta
bedtools getfasta -fi "${genome_fa}" -bed SNP_slop100.bed -fo SNP_slop100.fa -name
echo "got fasta"
date

### build REF and ALT fasta explicitly and compute SNP position within sequence
cp "${snp_bed}" SNP.bed
python3 <<'PY'
import sys

with open("SNP.bed") as f:
    fields = f.readline().strip().split()

if len(fields) < 6:
    raise ValueError("SNP.bed must contain at least 6 columns: chrom start end id ref alt")

chrom, start, end, snp_id, ref, alt = fields[:6]

start = int(start)
end   = int(end)
ref   = ref.upper()
alt   = alt.upper()

if "," in alt:
    raise ValueError(f"Multi-allelic variant not supported: ALT={alt}")

with open("SNP_slop100.bed") as f:
    wchrom, wstart, wend = f.readline().strip().split()[:3]

wstart = int(wstart)
wend   = int(wend)

variant_idx0 = start - wstart
variant_pos1 = variant_idx0 + 1

lines = [x.strip() for x in open("SNP_slop100.fa") if x.strip()]
seq = "".join(lines[1:]).upper()

if variant_idx0 < 0 or variant_idx0 >= len(seq):
    raise ValueError(f"Variant index out of range: variant_idx0={variant_idx0}, seq_len={len(seq)}")

fasta_ref_span = seq[variant_idx0:(variant_idx0 + len(ref))]
fasta_alt_span = seq[variant_idx0:(variant_idx0 + len(alt))]

if fasta_ref_span == ref:
    fasta_anchor_allele = "REF"
    fasta_anchor_len = len(ref)
    fasta_anchor_seq = fasta_ref_span
elif fasta_alt_span == alt:
    fasta_anchor_allele = "ALT"
    fasta_anchor_len = len(alt)
    fasta_anchor_seq = fasta_alt_span
    print(
        f"[WARN] FASTA matches BED ALT rather than BED REF; keeping BED REF/ALT labels: "
        f"bed_REF={ref}, bed_ALT={alt}, fasta_at_variant={fasta_anchor_seq}, "
        f"variant={chrom}:{start}-{end}",
        file=sys.stderr,
    )
else:
    raise ValueError(
        f"Neither BED allele matches FASTA at variant site: "
        f"bed_REF={ref}, fasta_REF_span={fasta_ref_span}, "
        f"bed_ALT={alt}, fasta_ALT_span={fasta_alt_span}, "
        f"variant={chrom}:{start}-{end}. "
        f"Check genome build, coordinates, and allele columns."
    )

left  = seq[:variant_idx0]
right = seq[(variant_idx0 + fasta_anchor_len):]

seq_ref = left + ref + right
seq_alt = left + alt + right

ref_start1 = variant_pos1
ref_end1   = variant_pos1 + len(ref) - 1

alt_start1 = variant_pos1
alt_end1   = variant_pos1 + len(alt) - 1

variant_type = (
    "SNV" if len(ref) == 1 and len(alt) == 1 else
    "deletion" if len(ref) > len(alt) else
    "insertion" if len(ref) < len(alt) else
    "MNV"
)

print(f"[DEBUG] variant={snp_id}")
print(f"[DEBUG] variant_type={variant_type}")
print(f"[DEBUG] window={wchrom}:{wstart}-{wend}")
print(f"[DEBUG] genomic={chrom}:{start}-{end}")
print(f"[DEBUG] original_seq_len={len(seq)}")
print(f"[DEBUG] REF_seq_len={len(seq_ref)}")
print(f"[DEBUG] ALT_seq_len={len(seq_alt)}")
print(f"[DEBUG] variant_idx0={variant_idx0}")
print(f"[DEBUG] variant_pos1={variant_pos1}")
print(f"[DEBUG] fasta_anchor_allele={fasta_anchor_allele}")
print(f"[DEBUG] fasta_anchor_at_variant={fasta_anchor_seq}")
print(f"[DEBUG] bed_REF={ref}")
print(f"[DEBUG] bed_ALT={alt}")

with open("SNP_slop100_REF.fa", "w") as out:
    out.write(f">{snp_id}|REF|variant_pos={variant_pos1}|variant_type={variant_type}\n")
    out.write(seq_ref + "\n")

with open("SNP_slop100_ALT.fa", "w") as out:
    out.write(f">{snp_id}|ALT|variant_pos={variant_pos1}|variant_type={variant_type}\n")
    out.write(seq_alt + "\n")

# Keep old file name for downstream compatibility.
with open("snp_position_in_seq.txt", "w") as out:
    out.write(str(variant_pos1) + "\n")

with open("variant_interval_REF.txt", "w") as out:
    out.write(f"{ref_start1}\t{ref_end1}\n")

with open("variant_interval_ALT.txt", "w") as out:
    out.write(f"{alt_start1}\t{alt_end1}\n")

with open("variant_info.tsv", "w") as out:
    out.write("variant_id\tchrom\tstart\tend\tref\talt\tvariant_type\tvariant_pos1\tref_start1\tref_end1\talt_start1\talt_end1\tref_seq_len\talt_seq_len\n")
    out.write(
        f"{snp_id}\t{chrom}\t{start}\t{end}\t{ref}\t{alt}\t{variant_type}\t"
        f"{variant_pos1}\t{ref_start1}\t{ref_end1}\t{alt_start1}\t{alt_end1}\t"
        f"{len(seq_ref)}\t{len(seq_alt)}\n"
    )
PY

variant_pos_seq=$(cat snp_position_in_seq.txt)
echo "Variant anchor position in sequence (1-based) = ${variant_pos_seq}"

### run FIMO on REF and ALT — no p-value threshold, capture all hits
for allele in REF ALT; do
    fimo --oc "results_${allele}" --thresh 1 --verbosity 1 "${memeF}" "SNP_slop100_${allele}.fa"

    ### strip comment lines and header
    awk 'BEGIN{FS=OFS="\t"} !/^#/ && NR>1 && $1!=""' \
        "results_${allele}/fimo.tsv" > "fimo_all_${allele}.tsv"

    ### keep only motifs whose hit window overlaps the SNP position
    ### FIMO columns: 1=motif_id 2=motif_alt_id 3=sequence_name
    ###               4=start 5=stop 6=strand 7=score 8=p-value 9=q-value 10=matched_sequence
    interval_file="variant_interval_${allele}.txt"
    var_start=$(cut -f1 "${interval_file}")
    var_end=$(cut -f2 "${interval_file}")

    awk -v snp="${snp_id}" -v allele="${allele}" \
        -v var_start="${var_start}" -v var_end="${var_end}" '
        BEGIN{FS=OFS="\t"}
        { if ($4 <= var_end && $5 >= var_start) print snp, allele, $0 }
    ' "fimo_all_${allele}.tsv" > "fimo_snp_overlap_${allele}.tsv"
done

echo "ran FIMO on REF and ALT"
date

### merge REF/ALT overlapping hits and classify differences
#source /home/shuo.shan-umw/miniconda3/etc/profile.d/conda.sh
#module load conda
#source activate fastQTL
/home/shuo.shan-umw/miniconda3/envs/fastQTL/bin/python3 <<'PY'
import pandas as pd
import numpy as np

def read_or_empty(fn):
    cols = [
        "snp", "allele", "motif_id", "motif_alt_id", "sequence_name",
        "start", "stop", "strand", "score", "pvalue", "qvalue", "matched_sequence"
    ]
    try:
        df = pd.read_csv(fn, sep="\t", header=None, names=cols)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=cols)
    return df

ref = read_or_empty("fimo_snp_overlap_REF.tsv")
alt = read_or_empty("fimo_snp_overlap_ALT.tsv")

def best_per_motif(df, suffix):
    if df.empty:
        return pd.DataFrame(columns=[
            "motif_id",
            f"start_{suffix}", f"stop_{suffix}", f"strand_{suffix}",
            f"score_{suffix}", f"pvalue_{suffix}", f"qvalue_{suffix}",
            f"matched_sequence_{suffix}"
        ])
    df = df.sort_values(["motif_id", "pvalue"], ascending=[True, True])
    df = df.groupby("motif_id", as_index=False).first()
    return df.rename(columns={
        "start":           f"start_{suffix}",
        "stop":            f"stop_{suffix}",
        "strand":          f"strand_{suffix}",
        "score":           f"score_{suffix}",
        "pvalue":          f"pvalue_{suffix}",
        "qvalue":          f"qvalue_{suffix}",
        "matched_sequence": f"matched_sequence_{suffix}"
    })[[
        "motif_id",
        f"start_{suffix}", f"stop_{suffix}", f"strand_{suffix}",
        f"score_{suffix}", f"pvalue_{suffix}", f"qvalue_{suffix}",
        f"matched_sequence_{suffix}"
    ]]

refb = best_per_motif(ref, "REF")
altb = best_per_motif(alt, "ALT")

merged = refb.merge(altb, on="motif_id", how="outer")

merged["present_REF"] = ~merged["score_REF"].isna()
merged["present_ALT"] = ~merged["score_ALT"].isna()

def classify(row):
    if     row["present_REF"] and not row["present_ALT"]: return "lost_in_ALT"
    elif not row["present_REF"] and row["present_ALT"]:   return "gained_in_ALT"
    elif   row["present_REF"] and     row["present_ALT"]:
        if   row["score_ALT"] > row["score_REF"]: return "stronger_in_ALT"
        elif row["score_ALT"] < row["score_REF"]: return "weaker_in_ALT"
        else:                                     return "unchanged"
    else:
        return "absent_both"

merged["neg_log10p_REF"]            = -np.log10(merged["pvalue_REF"].clip(lower=1e-300))
merged["neg_log10p_ALT"]            = -np.log10(merged["pvalue_ALT"].clip(lower=1e-300))
merged["delta_score_ALT_minus_REF"] = merged["score_ALT"]      - merged["score_REF"]
merged["delta_log10p_ALT_minus_REF"]= merged["neg_log10p_ALT"] - merged["neg_log10p_REF"]
merged["priority_abs_delta_log10p"] = merged["delta_log10p_ALT_minus_REF"].abs()
merged["motif_change_class"]        = merged.apply(classify, axis=1)

merged = merged.sort_values(
    ["motif_change_class", "priority_abs_delta_log10p"],
    ascending=[True, False],
    na_position="last"
)

merged.to_csv("fimo_REF_ALT_comparison.tsv", sep="\t", index=False)
PY

### move final output
[[ -s fimo_REF_ALT_comparison.tsv ]] || {
    echo "ERROR: merge failed; fimo_REF_ALT_comparison.tsv was not created" >&2
    exit 1
}
mv fimo_REF_ALT_comparison.tsv "${workingdir}/fimo_output_${snp_id}.txt"

echo "merged FIMO output written to ${workingdir}/fimo_output_${snp_id}.txt"
date

### clean up intermediates
rm -f SNP_slop100.bed SNP_slop100.fa \
      fimo_all_REF.tsv fimo_all_ALT.tsv \
      fimo_snp_overlap_REF.tsv fimo_snp_overlap_ALT.tsv
rm -rf "${workingdir}/results_REF" "${workingdir}/results_ALT"
echo "fimo already done for snp" > "${workingdir}/fimo_status.txt"


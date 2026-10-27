#!/bin/bash
# fetch chromBPNet variant score of all tested SNPs of a given gene, in all celltype x condition

SCRIPTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/script
OUTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats

# ---------- fetch input arguments ------------------- #
g=$1   # g=ITGA1
letter=${g:0:1}

mkdir -p "${OUTDIR}/${letter}/${g}"
cd "${OUTDIR}/${letter}/${g}"

echo "starting the script for ${g}"
date

out_file="${OUTDIR}/${letter}/${g}/chromBPNet_variant_scores.txt"
temp_out="${OUTDIR}/${letter}/${g}/chromBPNet_variant_scores.tmp.txt"

rm -f "$out_file" "$temp_out"

# store header in temp file
template_file="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_MEL_PBS_merged_alldonors/biasFactor0.5/variant_scores/summary.mean.variant_scores.tsv"

awk 'BEGIN { FS=OFS="\t" } NR == 1 { print "celltype", "condition", $0 }' \
    "$template_file" > "$temp_out"

# store SNP list
awk 'NR > 1' "${OUTDIR}/${letter}/${g}/genotype.txt" | cut -f3 > "${OUTDIR}/${letter}/${g}/snplist.txt"

# append matching chromBPNet variant scores to temp file
for celltype in MEL KRT FRB; do
    for condition in PBS IFNG; do

        if [ "$condition" = "IFNG" ]; then
            cond="IFN"
        elif [ "$condition" = "PBS" ]; then
            cond="PBS"
        fi

        folder="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_${celltype}_${cond}_merged_alldonors/biasFactor0.5/variant_scores"
        scoreFile="${folder}/summary.mean.variant_scores.tsv"

        awk -v ct="${celltype}" -v cond="${condition}" '
        BEGIN { FS=OFS="\t" }

        NR == FNR {
            snps[$1] = 1
            next
        }

        FNR > 1 {
            split($3, a, "_")
            if (a[1] in snps) {
                print ct, cond, $0
            }
        }
        ' "${OUTDIR}/${letter}/${g}/snplist.txt" "$scoreFile" >> "$temp_out"

        echo "$celltype" "$condition"
    done
done

# add SNP column to temp file and write final output file
awk '
BEGIN { FS=OFS="\t" }

NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "variant_id") variant_col = i
    }

    print "SNP", $0
    next
}

{
    snp = $variant_col
    sub(/_.*/, "", snp)

    print snp, $0
}
' "$temp_out" > "$out_file"

# remove intermediate file
rm -f "$temp_out"

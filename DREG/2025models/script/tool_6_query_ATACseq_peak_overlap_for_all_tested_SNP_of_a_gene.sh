#!/bin/bash
# query ATACseq peak overlap info of all tested SNP of a given gene

SCRIPTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/script
OUTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats

# ---------- fetch input arguments ------------------- #
g=$1   #g=ITGA1
letter=${g:0:1}
mkdir -p ${OUTDIR}/${letter}/${g}
cd ${OUTDIR}/${letter}/${g}
echo "starting the script for ${g}"; date

# turn genotype.txt into bed file
awk '
BEGIN { FS=OFS="\t" }
NR == 1 { next }
{
    start = $2 - 1
    end = $2
    print $1, start, end, $3
}
' genotype.txt > snp.bed

workingdir=${OUTDIR}/${letter}/${g}
SNPbed=${workingdir}/snp.bed
prefix=snp
out_file="${workingdir}/${prefix}_SNP_ATAC_overlap.txt"
rm -f ${out_file}
bash /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/02_overlap_atac.sh ${workingdir} ${prefix} ${SNPbed}


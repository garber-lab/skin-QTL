#!/bin/bash
# call Rscript to screen candidate chromatin-disrupting significant QTLs of a given gene from the big table compiled from tool7

SCRIPTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/script
OUTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats

# ---------- fetch input arguments ------------------- #
snp=$1
g=$2   #g=ITGA1
celltype=$3
condition=$4

#snp=rs12127787;g=GBP3;celltype=MEL;condition=IFNG
#snp=rs12127787;g=GBP3;celltype=KRT;condition=PBS
#snp=rs8141918;g=A4GALT;celltype=KRT;condition=PBS

letter=${g:0:1}
mkdir -p ${OUTDIR}/${letter}/${g}/fimo
mkdir -p ${OUTDIR}/${letter}/${g}/fimo/${snp}
workingdir=${OUTDIR}/${letter}/${g}/fimo/${snp}

#if [ -f "${workingdir}/fimo_summary_${snp}_${celltype}_${condition}.pdf" ]; then
#    echo "Output already exists: fimo_summary_${snp}_${celltype}_${condition}.pdf"
#    echo "Skipping ${g} : ${snp}"
#    exit 0
#fi

echo "starting the script for ${g} : ${snp}"; date
cd ${OUTDIR}/${letter}/${g}/fimo/${snp}

# ---- call fimo on variant ----
status_file="${workingdir}/fimo_status.txt"
fimo_output="${workingdir}/fimo_output_${snp}.txt"

grep -w "${snp}" "${OUTDIR}/${letter}/${g}/genotype.txt" \
    | awk 'BEGIN{OFS="\t"} {print $1,$2-1,$2,$3,$4,$5}' \
    > "${workingdir}/snp_bed6_format.bed"

bash /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/04_run_fimo.sh \
    "${workingdir}" "${snp}" "${workingdir}/snp_bed6_format.bed"

#if [[ ! -f "${status_file}" ]]; then
#   echo "fimo hasn't run for ${snp}, running now"; date
#
#    grep -w "${snp}" "${OUTDIR}/${letter}/${g}/genotype.txt" \
#        | awk 'BEGIN{OFS="\t"} {print $1,$2-1,$2,$3,$4,$5}' \
#        > "${workingdir}/snp_bed6_format.bed"
#
#    bash /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/04_run_fimo.sh \
#        "${workingdir}" "${snp}" "${workingdir}/snp_bed6_format.bed"
#
#else
#    if grep -q "fimo already done for snp" "${status_file}"; then
#        echo "fimo already done for snp, continue"
#        date
#
#    elif grep -q "fimo starting for snp" "${status_file}"; then
#        echo "fimo already running for snp ${snp}; waiting..."
#        date
#
#        until grep -q "fimo already done for snp" "${status_file}"; do
#            sleep 5
#        done
#
#        echo "fimo already done for snp, continue"
#        date
#    fi
#fi

# ---- call fimo result analysis script on variant ----
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
        Rscript /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/04_run_fimo_analyze.R  ${workingdir} ${snp} ${celltype} ${condition}
echo "done with fimo post-run analysis, plotting next."; date

# ---- make plots for important variants for fimo ----
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
        Rscript /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/04_compile_fimo_output_plots.R ${workingdir} ${snp} ${celltype} ${condition}

# ---- clean-up intermediate files ----
rm ${workingdir}/fimo_summary_${snp}_${celltype}_${condition}_objects.RData

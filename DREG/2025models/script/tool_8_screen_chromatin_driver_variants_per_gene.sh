#!/bin/bash
# call Rscript to screen candidate chromatin-disrupting significant QTLs of a given gene from the big table compiled from tool7

SCRIPTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/script
OUTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats

# ---------- fetch input arguments ------------------- #
g=$1   #g=ITGA1
letter=${g:0:1}
mkdir -p ${OUTDIR}/${letter}/${g}
cd ${OUTDIR}/${letter}/${g}
echo "starting the script for ${g}"; date

# ------------ call Rscript to plot --------------- #
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/tool_8_screen_chromatin_driver_variants_per_gene.R ${g} 

echo "written to ${OUTDIR}/${letter}/${g}/mechanistic_variant_table_detail.txt"
echo "written to ${OUTDIR}/${letter}/${g}/mechanistic_variant_list.txt"
echo "${g} has $(grep "rs" ${OUTDIR}/${letter}/${g}/mechanistic_variant_list.txt | wc -l | cut -d' ' -f1) candidate chromatin-disrupting QTL variants"

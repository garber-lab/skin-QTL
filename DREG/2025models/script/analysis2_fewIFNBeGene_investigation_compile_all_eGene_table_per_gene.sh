#!/bin/bash
# call Rscript to compile all modeling stats for all tested SNPs and significance annotation for each celltype x condition x QTLtype

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
        Rscript ${SCRIPTDIR}/analysis2_fewIFNBeGene_investigation_compile_all_eGene_table_per_gene.R ${g} ${OUTDIR}/${letter}/${g}


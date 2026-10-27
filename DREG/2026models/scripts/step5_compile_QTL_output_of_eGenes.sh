#!/bin/bash
# Compile, per context (celltype x condition x QTLtype), all tested variant-gene
# pairs belonging to significant eGenes (post-QC stats sliced to the fdr05 eGene
# set). Light I/O -- runs in a few minutes.
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8000]"
#BSUB -W 02:00
#BSUB -q short
#BSUB -J compileEgeneVars
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step5_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step5_%J_%I.err"

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models
SCRIPTDIR=${ROOT}/scripts

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step5_compile_QTL_output_of_eGenes.R ${ROOT}

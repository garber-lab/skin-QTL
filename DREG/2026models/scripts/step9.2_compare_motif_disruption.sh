#!/bin/bash
# STEP 9.2 -- call per-signal TF-motif disruption from allelic FIMO and compare
# co-signals per >=2-signal gene (different-TF vs same-TF; co-disrupted TF pairs).
# Run AFTER step9.1. Light.
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8000]"
#BSUB -W 01:00
#BSUB -q short
#BSUB -J fimo_compare
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step9.2_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step9.2_%J.err"

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models
SCRIPTDIR=${ROOT}/scripts
DELTA=2

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step9.2_compare_motif_disruption.R ${ROOT} ${DELTA}

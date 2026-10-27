#!/bin/bash
# STEP 8.2 -- pairwise cross-context signal matching: place each context's primary
# lead within the OTHER context's signal hierarchy (same-primary / rank-shift /
# distinct / novel-context-specific). Run AFTER step8.1. Writes:
#   results/conditional_analysis/pairwise/matched_<ctxA>_vs_<ctxB>.txt
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=48000]"
#BSUB -W 08:00
#BSUB -q long
#BSUB -J cond82_match
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step8.2_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step8.2_%J.err"

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models
SCRIPTDIR=${ROOT}/scripts
NEXPRPC=10; NGPC=3; LD_SAME=0.8

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step8.2_pairwise_signal_matching.R ${ROOT} ${NEXPRPC} ${NGPC} ${LD_SAME}

#!/bin/bash
# STEP 8.1 -- within-context stepwise (COJO-style) conditional signal discovery:
# for each context (celltype x condition) and each eGene, find the independent cis
# signals (primary, secondary, ...). Writes:
#   results/conditional_analysis/within_context/signals_<ct>_<cond>.txt
# Run FIRST (before step8.2). Loads the ~1.2GB genotype matrix -> high mem.
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=48000]"
#BSUB -W 08:00
#BSUB -q long
#BSUB -J cond81_signals
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step8.1_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step8.1_%J.err"

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models
SCRIPTDIR=${ROOT}/scripts
NEXPRPC=10; NGPC=3; MAXK=5; ALPHA=0.05   # covariates MUST match the main step1 run

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step8.1_within_context_stepwise_signals.R ${ROOT} ${NEXPRPC} ${NGPC} ${MAXK} ${ALPHA}

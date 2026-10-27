#!/bin/bash
# STEP 8.3 -- BH correction + finalize class + summary. Run AFTER step8.2. Writes:
#   results/conditional_analysis/pairwise/corrected/*, SUMMARY_per_comparison.txt,
#   INTERESTING_findings.txt
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8000]"
#BSUB -W 01:00
#BSUB -q short
#BSUB -J cond83_summary
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step8.3_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step8.3_%J.err"

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models
SCRIPTDIR=${ROOT}/scripts
FDR=0.05; MARG=0.05; LD_WEAK=0.5; PROX_BP=2000

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step8.3_correction_and_summary.R ${ROOT} ${FDR} ${MARG} ${LD_WEAK} ${PROX_BP}

#!/bin/bash
# Emit the master list of variants that are significant eQTLs OR reQTLs in ANY
# context (union over all eGenes), back-calculating per-variant significance from
# the pipeline's eigenMT(Sidak)+BH scheme. NO peak filter -- peak occupancy is a
# downstream ChromBPNet prioritization annotation, not an inclusion gate.
# Mirrors the sig-calling in chromBPNet_analysis/src/07_finemap_sig_eqtl.R.
# Writes: <ROOT>/results/list_of_variants_sig_eQTL_or_reQTL_of_union_eGenes.txt
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8000]"
#BSUB -W 01:00
#BSUB -q short
#BSUB -J sigQTLunionEgenes
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step5_sigQTL_union_eGenes_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step5_sigQTL_union_eGenes_%J.err"

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models
SCRIPTDIR=${ROOT}/scripts
ALPHA=0.05

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step5_get_sig_QTLs_of_union_eGenes.R ${ROOT} ${ALPHA}

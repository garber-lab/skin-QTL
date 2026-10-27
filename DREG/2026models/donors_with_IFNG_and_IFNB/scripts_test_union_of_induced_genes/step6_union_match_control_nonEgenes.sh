#!/bin/bash
# UNION variant of step6: build matched non-eGene controls for the UNION eGenes
# and compile their post-QC tested variants. Reads/writes
# <ROOT>/results_test_union_of_induced_genes (via the union .R, results_dir set there).
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=16000]"
#BSUB -W 03:00
#BSUB -q short
#BSUB -J matchControls_union
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/step6_union_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/step6_union_%J_%I.err"

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB
SCRIPTDIR=${ROOT}/scripts_test_union_of_induced_genes
BASELINE_SOURCE=CPM   # or VST

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step6_union_match_control_nonEgenes.R ${ROOT} ${BASELINE_SOURCE}

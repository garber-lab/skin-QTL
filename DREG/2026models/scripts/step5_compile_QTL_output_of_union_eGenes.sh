#!/bin/bash
# UNION-of-eGenes variant of step5 (genome-wide mapping): slice EVERY context to
# the GLOBAL union of eGenes (any gene that is an eGene in any celltype x any
# condition, eQTL), so every context has the cis profile of any eGene -- enabling
# colocalization ACROSS CELLTYPES and ACROSS CYTOKINES. Reuses the post-QC chunk
# tables (<ROOT>/<ct>/results_QC) -- unchanged. Writes, per context:
#   results/QTL_output_of_all_tested_variants_of_union_of_eGenes_<ct>_<cond>_eQTL.txt
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8000]"
#BSUB -W 02:00
#BSUB -q short
#BSUB -J compileUnionEgeneVars
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step5_union_eGenes_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/step5_union_eGenes_%J_%I.err"

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models
SCRIPTDIR=${ROOT}/scripts
QTL=eQTL

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step5_compile_QTL_output_of_union_eGenes.R ${ROOT} ${QTL}

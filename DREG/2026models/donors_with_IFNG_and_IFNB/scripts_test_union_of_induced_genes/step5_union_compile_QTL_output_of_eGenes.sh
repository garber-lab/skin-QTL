#!/bin/bash
# UNION variant of step5: compile all tested variant-gene pairs for the UNION
# eGenes (auto-discovered from results_test_union_of_induced_genes/*_gene_fdr05_genelist.txt).
# Reuses the SHARED post-QC chunk tables (<ROOT>/<ct>/results_QC) -- unchanged.
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8000]"
#BSUB -W 02:00
#BSUB -q short
#BSUB -J compileEgeneVars_union
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/step5_union_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/step5_union_%J_%I.err"

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB
SCRIPTDIR=${ROOT}/scripts_test_union_of_induced_genes

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step5_union_compile_QTL_output_of_eGenes.R ${ROOT}

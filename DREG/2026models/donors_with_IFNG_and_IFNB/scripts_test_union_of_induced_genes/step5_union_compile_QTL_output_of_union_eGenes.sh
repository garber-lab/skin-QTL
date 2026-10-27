#!/bin/bash
# UNION-eGenes variant of step5: for each celltype, compile all tested variant-gene
# pairs for the UNION of IFNG and IFNB eGenes (eQTL), in BOTH conditions -- so every
# union eGene has a full cis profile in both conditions and coloc can run on the
# IFNG-only / IFNB-only categories, not just shared. Reuses the SHARED post-QC chunk
# tables (<ROOT>/<ct>/results_QC) -- unchanged. Writes:
#   results_test_union_of_induced_genes/
#     QTL_output_of_all_tested_variants_of_UNION_eGenes_<ct>_<cond>_eQTL.txt
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8000]"
#BSUB -W 02:00
#BSUB -q short
#BSUB -J compileUnionEgeneVars
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/step5_union_eGenes_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/step5_union_eGenes_%J_%I.err"

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB
SCRIPTDIR=${ROOT}/scripts_test_union_of_induced_genes

# condition pair + QTL type to union over (defaults inside the R script)
COND_A=IFNG
COND_B=IFNB
QTL=eQTL

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step5_union_compile_QTL_output_of_union_eGenes.R ${ROOT} ${COND_A} ${COND_B} ${QTL}

#!/bin/bash
# Build a matched non-eGene control set (non-eGenes matched to each eGene on
# baseline PBS CPM + DESeq2 log2FC, with RELATIVE calipers on lead-SNP min-genotype
# count and n_snps_used), compile the control set's post-QC tested variants, and
# write per-context matching + diagnostics tables. Tune the calipers/q-buffer at
# the top of the .R (POOL_QMIN / TOL_POWER / TOL_NSNPS); the diagnostics section
# recommends new values if any caliper is too tight.
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=16000]"
#BSUB -W 03:00
#BSUB -q short
#BSUB -J matchControls
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNB_and_TNF/logs/step6_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNB_and_TNF/logs/step6_%J_%I.err"

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNB_and_TNF
SCRIPTDIR=${ROOT}/scripts
BASELINE_SOURCE=CPM   # or VST

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step6_match_control_nonEgenes.R ${ROOT} ${BASELINE_SOURCE}

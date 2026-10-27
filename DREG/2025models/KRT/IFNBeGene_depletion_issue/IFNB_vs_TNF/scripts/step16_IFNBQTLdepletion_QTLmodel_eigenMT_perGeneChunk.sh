#!/bin/bash
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -W 1:00
#BSUB -q short
#BSUB -J QTL.KRT[1-80]%80
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/KRT/IFNBeGene_depletion_issue/IFNB_vs_TNF/logs/eigenMT_QTL_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/KRT/IFNBeGene_depletion_issue/IFNB_vs_TNF/logs/eigenMT_QTL_%J_%I.err"

sleep $((RANDOM % 20)) # reduces heavy I/O thundering-herd effect in parallel jobs
ct=KRT
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}/IFNBeGene_depletion_issue/IFNB_vs_TNF
SCRIPTDIR=${DIR}/scripts
mkdir -p ${DIR}/eigenMT
mkdir -p ${DIR}/eigenMT/TNF/eQTL
mkdir -p ${DIR}/eigenMT/IFNB/eQTL

# Convert LSB_JOBINDEX → chunk_id (e.g. 1-120 → 001–120)
chunk_id=$(printf "%03d" $((LSB_JOBINDEX)))

## ---------------------------------------------------------- #
# Submit chunk to eigenMT Rscript, it calculates FWER adjusted p.val for gene:SNP pairs that passed QC
echo "Performing eigenMT for ${chunk_id}";date
cd ${DIR}

# Tell R to use my library inside the container:
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step16_IFNBQTLdepletion_QTLmodel_eigenMT_perGeneChunk.R ${ct} TNF IFNB eQTL ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTDIR}/step16_IFNBQTLdepletion_QTLmodel_eigenMT_perGeneChunk.R ${ct} IFNB TNF eQTL ${chunk_id}

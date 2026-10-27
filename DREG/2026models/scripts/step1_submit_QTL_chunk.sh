#!/bin/bash
#BSUB -n 5
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4000]"
#BSUB -W 72:00
#BSUB -q long
#BSUB -J QTL[250,262,270,292]
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/qtl2_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/qtl2_%J_%I.err"
### QTL[1-314]%200 usually -R "rusage[mem=2000]"

sleep $((RANDOM % 20)) # reduces heavy I/O thundering-herd effect in parallel jobs

DATADIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/data
SCRIPTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/scripts
mkdir -p /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/FRB/results
mkdir -p /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/MEL/results
mkdir -p /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/KRT/results

# ---- Select chunk file safely ---- #
# e.g. chunk=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/data/chunks/pairs_chunk_001.tsv; chunk_base=pairs_chunk_001; chunk_id=001
chunk=$(ls ${DATADIR}/chunks/pairs_chunk_*.tsv | sort -V | sed -n ${LSB_JOBINDEX}p)
chunk_base=$(basename "$chunk" .tsv)
chunk_id=$(echo "$chunk_base" | sed 's/pairs_chunk_//')

# -------- Set-up Singularity and Dependencies for Rscript ------ #
echo "Fitting QTL model for ${chunk}";date

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
# Per-cell-type covariate settings chosen from the optimizer:
#   step1_run_QTL_chunk.R <CELLTYPE> <nExprPC (condition-removed)> <nGPC> <chunk_id>
#singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif Rscript ${SCRIPTDIR}/step1_run_QTL_chunk.R FRB 10 3 ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif Rscript ${SCRIPTDIR}/step1_run_QTL_chunk.R KRT 10 3 ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif Rscript ${SCRIPTDIR}/step1_run_QTL_chunk.R MEL 10 3 ${chunk_id}

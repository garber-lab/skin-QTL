#!/bin/bash
#BSUB -n 5
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -W 24:00
#BSUB -q long
#BSUB -J QTL[1-94]%94
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/qtl_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/qtl_%J_%I.err"
### Set the array upper bound to the chunk count reported by step0 (Num gene chunks).

sleep $((RANDOM % 20)) # reduces heavy I/O thundering-herd effect in parallel jobs

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB   # this pipeline's outputs + its own chunks
DATADIR=${ROOT}/data
SCRIPTDIR=${ROOT}/scripts
mkdir -p ${ROOT}/FRB/results
mkdir -p ${ROOT}/MEL/results
mkdir -p ${ROOT}/KRT/results

# ---- Select chunk file safely ---- #
# e.g. chunk=${ROOT}/data/chunks/pairs_chunk_001.tsv; chunk_base=pairs_chunk_001; chunk_id=001
chunk=$(ls ${DATADIR}/chunks/pairs_chunk_*.tsv | sort -V | sed -n ${LSB_JOBINDEX}p)
chunk_base=$(basename "$chunk" .tsv)
chunk_id=$(echo "$chunk_base" | sed 's/pairs_chunk_//')

# -------- Set-up Singularity and Dependencies for Rscript ------ #
echo "Fitting QTL model (matched-donor: IFNG & IFNB) for ${chunk}";date

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
# Per-cell-type covariate settings chosen from the optimizer:
#   step1_run_QTL_chunk.R <CELLTYPE> <nExprPC (condition-removed)> <nGPC> <chunk_id>
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif Rscript ${SCRIPTDIR}/step1_run_QTL_chunk.R FRB 10 3 ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif Rscript ${SCRIPTDIR}/step1_run_QTL_chunk.R KRT 10 3 ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif Rscript ${SCRIPTDIR}/step1_run_QTL_chunk.R MEL 10 3 ${chunk_id}

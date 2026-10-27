#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -W 2:00
#BSUB -q short
#BSUB -J QC[1-94]
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/qc_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/qc_%J_%I.err"
### Set the array upper bound to the chunk count reported by step0 (Num gene chunks).

sleep $((RANDOM % 10)) # reduces heavy I/O thundering-herd effect in parallel jobs
ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB   # this pipeline's outputs + its own chunks
SCRIPTDIR=${ROOT}/scripts
mkdir -p ${ROOT}/FRB/results_QC
mkdir -p ${ROOT}/MEL/results_QC
mkdir -p ${ROOT}/KRT/results_QC

# Select chunk file (this pipeline's own chunks, built by step0)
# map 1..N -> 001..N
idx0=$((LSB_JOBINDEX))
chunk_id=$(printf "%03d" "${idx0}")
chunk="${ROOT}/data/chunks/pairs_chunk_${chunk_id}.tsv"

echo "JOBINDEX=${LSB_JOBINDEX} -> chunk_id=${chunk_id} -> ${chunk}"

# Call Rscript to run QC
echo "Performing QC (matched-donor: IFNG & IFNB) for ${chunk}";date

# Tell R to use my library inside the container:
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
bash -lc "
for ct in KRT MEL FRB; do
	for cond in PBS IFNG IFNB; do
		Rscript ${SCRIPTDIR}/step2_QTL_result_QC.R \$ct \$cond eQTL  ${chunk_id}
	done

	for cond in IFNG IFNB; do
		Rscript ${SCRIPTDIR}/step2_QTL_result_QC.R \$ct \$cond reQTL ${chunk_id}
	done
done
"

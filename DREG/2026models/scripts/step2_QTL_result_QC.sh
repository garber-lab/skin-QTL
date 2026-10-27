#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -W 2:00
#BSUB -q short
#BSUB -J QC[1-314]%200
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/qc_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/qc_%J_%I.err"

sleep $((RANDOM % 10)) # reduces heavy I/O thundering-herd effect in parallel jobs
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models
mkdir -p ${DIR}/FRB/results_QC
mkdir -p ${DIR}/MEL/results_QC
mkdir -p ${DIR}/KRT/results_QC

# Select chunk file
# e.g. chunk id is exactly the array index (0..188)
# e.g. map 1..189 -> 001..189
idx0=$((LSB_JOBINDEX))
chunk_id=$(printf "%03d" "${idx0}")
chunk="${DIR}/data/chunks/pairs_chunk_${chunk_id}.tsv"
chunk_base="pairs_chunk_${chunk_id}"

echo "JOBINDEX=${LSB_JOBINDEX} -> chunk_id=${chunk_id} -> ${chunk}"

# Call Rscript to run QC
echo "Performing QC for ${chunk}";date

# Tell R to use my library inside the container:
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
bash -lc "
for ct in KRT MEL FRB; do
	for cond in PBS IFNG IFNB TNF; do
		Rscript /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/scripts/step2_QTL_result_QC.R \$ct \$cond eQTL  ${chunk_id}
	done
	
	for cond in IFNG IFNB TNF; do
		Rscript /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/scripts/step2_QTL_result_QC.R \$ct \$cond reQTL ${chunk_id}
	done
done
"

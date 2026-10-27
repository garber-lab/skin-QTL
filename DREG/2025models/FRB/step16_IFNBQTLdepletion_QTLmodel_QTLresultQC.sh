#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -W 2:00
#BSUB -q short
#BSUB -J QTL.FRB[1-119]%119
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/IFNBeGene_depletion_issue/logs/qc_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/IFNBeGene_depletion_issue/logs/qc_%J_%I.err"

sleep $((RANDOM % 10)) # reduces heavy I/O thundering-herd effect in parallel jobs
ct=FRB
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}/IFNBeGene_depletion_issue
mkdir -p ${DIR}/results_QC

# Select chunk file
# chunk id is exactly the array index (0..188)
# map 1..189 -> 001..189
idx0=$((LSB_JOBINDEX))
chunk_id=$(printf "%03d" "${idx0}")
chunk="${DIR}/chunks/pairs_chunk_${chunk_id}.tsv"
chunk_base="pairs_chunk_${chunk_id}"

echo "JOBINDEX=${LSB_JOBINDEX} -> chunk_id=${chunk_id} -> ${chunk}"

# Call Rscript to run QC
echo "Performing QC for ${chunk}";date

# Tell R to use my library inside the container:
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
bash -lc "
for cond in PBS IFNG IFNB TNF; do
  Rscript /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}/step16_IFNBQTLdepletion_QTLmodel_QTLresultQC.R ${ct} \$cond eQTL  ${chunk_id}
done

for cond in IFNG IFNB TNF; do
  Rscript /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}/step16_IFNBQTLdepletion_QTLmodel_QTLresultQC.R ${ct} \$cond reQTL ${chunk_id}
done
"

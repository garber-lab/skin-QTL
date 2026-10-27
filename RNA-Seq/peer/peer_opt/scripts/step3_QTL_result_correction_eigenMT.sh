#!/bin/bash
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -W 2:00
#BSUB -q short
#BSUB -J eigenMT[1-314]%200
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/eigenMT_QTL_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/eigenMT_QTL_%J_%I.err"
# takes <5 minutes for run 1 chunk

sleep $((RANDOM % 20)) # reduces heavy I/O thundering-herd effect in parallel jobs
SCRIPTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/scripts

for ct in KRT MEL FRB; do
	for cond in PBS IFNG IFNB TNF; do
		DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/${ct}
		mkdir -p ${DIR}/eigenMT
		mkdir -p ${DIR}/eigenMT/${cond}/eQTL
		mkdir -p ${DIR}/eigenMT/${cond}/reQTL
		rm -rf ${DIR}/eigenMT/PBS/reQTL
	done
done

# Convert LSB_JOBINDEX → chunk_id (e.g. 1-120 → 001–120)
chunk_id=$(printf "%03d" $((LSB_JOBINDEX)))

## ---------------------------------------------------------- #
# Submit chunk to eigenMT Rscript, it calculates FWER adjusted p.val for gene:SNP pairs that passed QC
echo "Performing eigenMT for ${chunk_id}";date
cd ${DIR}

# Tell R to use my library inside the container:
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
bash -lc "
for ct in KRT MEL FRB; do
        for cond in PBS IFNG IFNB TNF; do
                Rscript /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/scripts/step3_QTL_result_correction_eigenMT.R \$ct \$cond eQTL  ${chunk_id}
        done

        for cond in IFNG IFNB TNF; do
                Rscript /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/scripts/step3_QTL_result_correction_eigenMT.R \$ct \$cond reQTL ${chunk_id}
        done
done
"

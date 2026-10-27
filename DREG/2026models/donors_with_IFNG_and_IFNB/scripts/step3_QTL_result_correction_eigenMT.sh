#!/bin/bash
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4000]"
#BSUB -W 2:00
#BSUB -q short
#BSUB -J eigenMT[1-94]
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/eigenMT_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/eigenMT_%J_%I.err"
# takes <5 minutes to run 1 chunk
### Set the array upper bound to the chunk count reported by step0 (Num gene chunks).

sleep $((RANDOM % 20)) # reduces heavy I/O thundering-herd effect in parallel jobs
ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB   # this pipeline's outputs
SCRIPTDIR=${ROOT}/scripts

for ct in KRT MEL FRB; do
	for cond in PBS IFNG IFNB; do
		DIR=${ROOT}/${ct}
		mkdir -p ${DIR}/eigenMT
		mkdir -p ${DIR}/eigenMT/${cond}/eQTL
		mkdir -p ${DIR}/eigenMT/${cond}/reQTL
	done
done

# Convert LSB_JOBINDEX → chunk_id (e.g. 1-314 → 001–314)
chunk_id=$(printf "%03d" $((LSB_JOBINDEX)))

## ---------------------------------------------------------- #
# Submit chunk to eigenMT Rscript, it calculates FWER adjusted p.val for gene:SNP pairs that passed QC
echo "Performing eigenMT (matched-donor: IFNG & IFNB) for ${chunk_id}";date

# Tell R to use my library inside the container:
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
bash -lc "
for ct in KRT MEL FRB; do
        for cond in PBS IFNG IFNB; do
                Rscript ${SCRIPTDIR}/step3_QTL_result_correction_eigenMT.R \$ct \$cond eQTL  ${chunk_id}
        done

        for cond in IFNG IFNB; do
                Rscript ${SCRIPTDIR}/step3_QTL_result_correction_eigenMT.R \$ct \$cond reQTL ${chunk_id}
        done
done
"

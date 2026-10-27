#!/bin/bash
# perform BH correction on locally adjusted pvalues across all genes
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -W 02:00
#BSUB -q short
#BSUB -J BH
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/BH_QTL_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/BH_QTL_%J_%I.err"

SCRIPTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/scripts
mkdir -p ${DIR}/eigenMT/results
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

# concatenate all eigenMT result tables 
for ct in MEL KRT FRB; do
for cond in PBS IFNG IFNB TNF; do

	DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/${ct}
	mkdir -p ${DIR}/eigenMT/results

	# --- eQTL ---
	echo "concatenating eigenMT results for ${ct} ${cond} eQTLs"
	cd ${DIR}/eigenMT/${cond}/eQTL
	head -1 ${ct}_${cond}_eQTL_001.eigenMT.tsv > ${DIR}/eigenMT/results/${ct}_${cond}_eQTL.eigenMT.txt
	for f in ${ct}_${cond}_eQTL*.tsv;do
		awk 'NR>1' ${f} >> ${DIR}/eigenMT/results/${ct}_${cond}_eQTL.eigenMT.txt
	done
	
	echo "doing BH correction for ${cond} eQTL across genes"
	cd ${DIR}/eigenMT/results
	singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
        Rscript ${SCRIPTDIR}/step4_QTL_result_correction_BH.R ${ct}_${cond}_eQTL.eigenMT.txt

	# --- reQTL ---
	if [ "${cond}" != "PBS" ]; then
		echo "concatenating eigenMT results for ${ct} ${cond} reQTLs"
		cd ${DIR}/eigenMT/${cond}/reQTL
		head -1 ${ct}_${cond}_reQTL_001.eigenMT.tsv > ${DIR}/eigenMT/results/${ct}_${cond}_reQTL.eigenMT.txt	
		for f in ${ct}_${cond}_reQTL*.tsv;do
			awk 'NR>1' ${f} >> ${DIR}/eigenMT/results/${ct}_${cond}_reQTL.eigenMT.txt
		done

		echo "doing BH correction for ${ct} ${cond} reQTL across genes"
		cd ${DIR}/eigenMT/results
		singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
        	Rscript ${SCRIPTDIR}/step4_QTL_result_correction_BH.R ${ct}_${cond}_reQTL.eigenMT.txt
	fi
done
done

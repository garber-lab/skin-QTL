#!/bin/bash
# perform BH correction on locally adjusted pvalues across all genes
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -W 01:00
#BSUB -q short
#BSUB -J BH
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/IFNBeGene_depletion_issue/logs/BH_QTL_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/IFNBeGene_depletion_issue/logs/BH_QTL_%J_%I.err"

ct=FRB
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}/IFNBeGene_depletion_issue
SCRIPTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
mkdir -p ${DIR}/eigenMT/results

# concatenate all eigenMT result tables 
for cond in IFNG IFNB;do
	echo "concatenating eigenMT results for ${cond} eQTLs"
	cd ${DIR}/eigenMT/${cond}/eQTL
	head -1 ${ct}_${cond}_eQTL_001.eigenMT.tsv > ${DIR}/eigenMT/results/${ct}_${cond}_eQTL.eigenMT.txt
	for f in ${ct}_${cond}_eQTL*.tsv;do
		awk 'NR>1' ${f} >> ${DIR}/eigenMT/results/${ct}_${cond}_eQTL.eigenMT.txt
	done
done

# perform BH
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

for cond in IFNG IFNB;do
	echo "doing BH correction for ${cond} eQTL across genes"
	cd ${DIR}/eigenMT/results
	singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
        Rscript ${SCRIPTDIR}/step16_IFNBQTLdepletion_QTLmodel_BHcorrection.R ${ct}_${cond}_eQTL.eigenMT.txt
done

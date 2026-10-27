#!/bin/bash
# UNION-FAMILY BH correction for the IFNG-vs-IFNB experiment.
# Redefines the IFNG & IFNB BH family as the UNION of the two induced sets and
# writes eGene calls into <ROOT>/results_test_union_of_induced_genes (the gated
# results dir is left untouched).
#
# Reuses the SAME concatenated eigenMT tables as the gated run
# (<ROOT>/<ct>/eigenMT/results/*.eigenMT.txt): these are family-independent, so
# steps 1-3 do NOT need to be re-run. If a concatenated table is missing it is
# rebuilt here (idempotent).
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -W 02:00
#BSUB -q short
#BSUB -J BH_union
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/BH_union_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB/logs/BH_union_%J_%I.err"

ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNG_and_IFNB
SCRIPTDIR=${ROOT}/scripts_test_union_of_induced_genes   # union scripts live here on HPC
RESULTS=${ROOT}/results_test_union_of_induced_genes
mkdir -p ${RESULTS}
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

# IFNG vs IFNB experiment only (PBS/TNF unchanged in the gated run)
for ct in MEL KRT FRB; do
for cond in IFNG IFNB; do

	DIR=${ROOT}/${ct}
	mkdir -p ${DIR}/eigenMT/results

	for qtl in eQTL reQTL; do
		AGG=${DIR}/eigenMT/results/${ct}_${cond}_${qtl}.eigenMT.txt

		# rebuild the concatenated eigenMT table only if missing (idempotent)
		if [ ! -s "${AGG}" ]; then
			echo "concatenating eigenMT results for ${ct} ${cond} ${qtl}"
			cd ${DIR}/eigenMT/${cond}/${qtl}
			head -1 ${ct}_${cond}_${qtl}_001.eigenMT.tsv > ${AGG}
			for f in ${ct}_${cond}_${qtl}*.tsv; do awk 'NR>1' ${f} >> ${AGG}; done
		else
			echo "reusing existing eigenMT table ${AGG}"
		fi

		echo "UNION BH correction for ${ct} ${cond} ${qtl}"
		cd ${DIR}/eigenMT/results
		singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
			Rscript ${SCRIPTDIR}/step4_union_BH.R ${ct}_${cond}_${qtl}.eigenMT.txt
	done
done
done

echo "done. UNION eGene calls in ${RESULTS}"

#!/bin/bash
# use the data from data/all_QTL_stats, compile useful information table for any gene that has a sig QTL in any condition

# set-up
mainDir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/fewIFNBeGene_investigation
mkdir -p ${mainDir}/logs

# compile list of genes
outFile_long=${mainDir}/eGene_list_long.txt
outFile_collapsed=${mainDir}/eGene_list.txt

rm -f ${outFile_long}
rm -f ${outFile_collapsed}

for ct in MEL KRT FRB;do
	cd /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}/eigenMT/results

	for f in *gene_fdr05_genelist.txt;do
		[[ "${f}" == *_PBS_eQTL_gene_fdr05_genelist.txt ]] && continue

		prefix=$(echo ${f} | cut -d'_' -f1-3)
		awk -v prefix=${prefix} '{OFS="\t"}{print $0, prefix}' ${f} >> ${outFile_long}
		echo ${f}
	done
	echo ${ct}
done

awk 'BEGIN{OFS="\t"} {a[$1]=($1 in a ? a[$1] "," $2 : $2)} END{for (g in a) print g, a[g]}' "${outFile_long}" | sort -k1 > ${outFile_collapsed}

# ---
scriptDir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/script
inDir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats
ls -d ${inDir}/*/* > ${inDir}/gene_dirs.txt

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
while read line;do

	g=$(echo ${line} | cut -d' ' -f1)
	letter=${g:0:1}

	g=GBP3;letter=${g:0:1}
	singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${scriptDir}/analysis2_fewIFNBeGene_investigation_compile_all_eGene_table_per_gene.R ${g} ${inDir}/${letter}/${g}



	#cat ${inDir}/${letter}/${g}/leads.txt | awk 'NR==1 || FNR>1' > ${inDir}/${letter}/${g}/master_leads.txt
	## enrichment test
	#singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
        #Rscript ${scriptDir}/analysis2_fewIFNBeGene_investigation_test_lead_cbp_enrichment.R ${inDir}/${letter}/${g}/master_leads.txt ${inDir}/gene_dirs.txt 

done < ${outFile_collapsed}


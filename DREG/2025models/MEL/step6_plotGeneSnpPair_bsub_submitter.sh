#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=1000]
#BSUB -q long
#BSUB -W 08:00
#BSUB -J job_submitter
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/job_submitter_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/job_submitter_%J.err"
# batch submit jobs to plot gene:snp pairs

# reduces heavy I/O thundering-herd effect in parallel jobs
sleep $((RANDOM % 20))

# set-up
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL
f=MEL_PBS_eQTL_gene_fdr05_table.txt
basename=$(echo ${f} | sed 's/\_table\.txt//g')
cd ${DIR}/plots

# create batch job commands
echo "creating commands for ${f}"
awk 'NR>1' ${DIR}/eigenMT/results/${f} | cut -f4,8 | tr "\t" " " > ${DIR}/plots/pairs_for_${f}
awk '{print "bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/step6_plotGeneSnpPair.sh",$0}' ${DIR}/plots/pairs_for_${f} > ${DIR}/plots/commands_for_${f}
#head -20 ${DIR}/plots/commands_for_${f} > ${DIR}/plots/top20.commands_for_${f}

# first sort gene name's first characters
echo "submitting $(wc -l ${DIR}/plots/commands_for_${f} | cut -d' ' -f1) jobs"; date
while read c; do
	echo ${c} | bsub -W 0:30 -J plot -n 1 -R "span[hosts=1]" -R rusage[mem=500] -q short -o "${DIR}/logs/plot_%J_${f}.out" -e "${DIR}/logs/plot_%J_${f}.err"
done < ${DIR}/plots/commands_for_${f}

# monitor jobs' progress
while [[ $(bjobs | grep 'plot' | wc -l) != 0 ]] ; do echo $(bjobs | grep 'plot' | wc -l) "jobs remaining"; sleep 5; done

# compile all PDF into one big PDF
cd ${DIR}/plots/temp_output
mkdir -p ${basename}
mv *.pdf ${basename}

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${DIR}/step6_combine_all_plots.R ${DIR}/plots/temp_output/${basename}

rm -r ${basename}

#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=1000]
#BSUB -q long
#BSUB -W 08:00
#BSUB -J job_submitter
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/fewIFNBeGene_investigation/logs/plotting_job_submitter_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/fewIFNBeGene_investigation/logs/plotting_job_submitter_%J.err"
# batch submit jobs to plot gene:snp pairs


mainDir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/fewIFNBeGene_investigation
mkdir -p ${mainDir}/logs
mkdir -p ${mainDir}/plots

# set-up
ct=MEL # FRB, KRT, MEL
inFile=${mainDir}/temp_${ct}_gene_snp_pairs_for_plotting.txt
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
basename=$(echo ${inFile} | sed 's/.*\///g' | sed 's/\.txt//g')
cd ${DIR}/plots
rm -f ${DIR}/plots/temp_output/*

# create batch job commands
echo "creating commands for ${inFile}"
awk -v celltype=${ct} '{print "bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/"celltype"/step6_plotGeneSnpPair.sh",$0}' ${inFile} > ${DIR}/plots/commands_for_${basename}

# first sort gene name's first characters
echo "submitting $(wc -l ${DIR}/plots/commands_for_${basename} | cut -d' ' -f1) jobs"; date
while read c; do
        echo ${c} | bsub -W 0:30 -J plot -n 1 -R "span[hosts=1]" -R rusage[mem=500] -q short -o "${mainDir}/logs/plot_%J_${basename}.out" -e "${mainDir}/logs/plot_%J_${basename}.err"
done < ${DIR}/plots/commands_for_${basename}

# monitor jobs' progress
while [[ $(bjobs | grep 'plot' | wc -l) != 0 ]] ; do echo $(bjobs | grep 'plot' | wc -l) "jobs remaining"; sleep 5; done

# move plots to output folder
cd ${DIR}/plots/temp_output
for f in *.pdf;do
	mv ${f} ${mainDir}/plots/${ct}_${f}
done



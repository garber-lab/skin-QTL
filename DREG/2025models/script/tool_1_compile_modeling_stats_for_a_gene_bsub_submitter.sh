#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -W 24:00
#BSUB -q long
#BSUB -J compile_stats
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats/bsub_job_submitter.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats/bsub_job_submitter.err"

# bsub job submitter
out_dir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats
mkdir -p ${out_dir}/logs

# list of genes
gene_file=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_genes.bed
gene_list=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_genes.txt
cut -f4 ${gene_file} | awk 'NR>1' | sort -u > ${gene_list}

# turn bsub job into array 
job_array_file=${out_dir}/job_array.txt
rm -f ${job_array_file}
awk '{print "bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/script/tool_compile_modeling_stats_for_a_gene.sh", $0}' ${gene_list} > ${job_array_file}

# collapse jobs into 50 commands per command
collapsed_job_array_file=${job_array_file%.txt}.collapsed50.txt

awk '{
  printf "%s%s", (NR % 50 == 1 ? "" : "; "), $0
  if (NR % 50 == 0) printf "\n"
}
END {
  if (NR % 50 != 0) printf "\n"
}' ${job_array_file} > ${collapsed_job_array_file}


# submit as bsub jobs
while read c; do
	echo ${c} | bsub -J compile_stats -W 8:00 -n 1 -R "span[hosts=1]" -R rusage[mem=1080] -q long -o "${out_dir}/logs/%J%I.out" -e "${out_dir}/logs/%J%I.err"
done < ${collapsed_job_array_file}


# clean-up
rm ${gene_list}
rm ${collapsed_job_array_file}
rm ${job_array_file}


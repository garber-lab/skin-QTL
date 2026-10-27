#!/bin/bash
# workflow to run FIMO on all shortlisted chromatin-disrupting variants for all genes (output from tool8)

SCRIPTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/script
OUTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats

# ---- fetch all variants and their genes ----
# list of genes
gene_file=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_genes.bed
gene_list=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_genes.txt
cut -f4 ${gene_file} | awk 'NR>1' | sort -u > ${gene_list}

variant_gene_dict_detail=${OUTDIR}/candidate_variant_gene_dict_detail.txt
rm -f ${variant_gene_dict_detail}
while IFS= read -r g; do

        letter=${g:0:1}
	variant_file=${OUTDIR}/${letter}/${g}/mechanistic_variant_table_detail.txt

        if [[ -s "${variant_file}" ]] && grep -q "rs" "${variant_file}"; then
                awk 'BEGIN {FS=OFS="\t"} /rs/ {print $1, $6, $14, $15}' "${variant_file}" >> "${variant_gene_dict_detail}"
        fi

done < ${gene_list}

# ---- generate BSUB job commands ----
# turn bsub job into array
job_array_file=${OUTDIR}/job_array_tool9.txt
rm -f ${job_array_file}
awk '{print "bash /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/script/tool_9_run_FIMO_on_candidate_variants_per_gene.sh", $1, $2, $3, $4}' ${variant_gene_dict_detail} > ${job_array_file} 

# collapse jobs into 100 commands per command
collapsed_job_array_file=${job_array_file%.txt}.collapsed100.txt
awk '{
  printf "%s%s", (NR % 100 == 1 ? "" : "; "), $0
  if (NR % 100 == 0) printf "\n"
}
END {
  if (NR % 100 != 0) printf "\n"
}' ${job_array_file} > ${collapsed_job_array_file}

# submit as bsub jobs
while read c; do
        echo ${c} | bsub -J compile_stats -W 72:00 -n 1 -R "span[hosts=1]" -R rusage[mem=4000] -q long -o "${OUTDIR}/logs/%J%I.out" -e "${OUTDIR}/logs/%J%I.err"
done < ${collapsed_job_array_file}



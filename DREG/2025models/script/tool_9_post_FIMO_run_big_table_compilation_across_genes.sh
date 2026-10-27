#!/bin/bash
# workflow to run FIMO on all shortlisted chromatin-disrupting variants for all genes (output from tool8)

SCRIPTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/script
OUTDIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats

# ---- fetch all variants and their genes ----
gene_file=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_genes.bed
gene_list=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_genes.txt
cut -f4 ${gene_file} | awk 'NR>1' | sort -u > ${gene_list}

fimo_output_masterfile=${OUTDIR}/fimo_output_masterfile.txt
rm -f ${fimo_output_masterfile}
header_written=0

while IFS= read -r g; do

    letter=${g:0:1}
    fimo_dir=${OUTDIR}/${letter}/${g}/fimo

    # skip genes with no fimo directory
    [[ ! -d "${fimo_dir}" ]] && echo "[SKIP] No fimo dir for ${g}" && continue

    # per-gene output (optional but useful)
    gene_outfile=${OUTDIR}/${letter}/${g}/fimo_output_masterfile_${g}.txt
    gene_header_written=0

    for variant_dir in ${fimo_dir}/rs*/; do
        [[ ! -d "${variant_dir}" ]] && continue
        variant=$(basename "${variant_dir}")

        for summary_file in ${variant_dir}/fimo_summary_${variant}_*_table.txt; do
            [[ ! -f "${summary_file}" ]] && continue

            if [[ ${header_written} -eq 0 ]]; then
                # write header with two new leading columns
                echo -e "gene\tvariant\t$(head -1 ${summary_file})" > ${fimo_output_masterfile}
                header_written=1
            fi

            if [[ ${gene_header_written} -eq 0 ]]; then
                echo -e "gene\tvariant\t$(head -1 ${summary_file})" > ${gene_outfile}
                gene_header_written=1
            fi

            # append data rows (skip header) with gene and variant prepended
            awk -v g="${g}" -v v="${variant}" \
                'NR>1 { print g "\t" v "\t" $0 }' \
                "${summary_file}" >> ${fimo_output_masterfile}

            awk -v g="${g}" -v v="${variant}" \
                'NR>1 { print g "\t" v "\t" $0 }' \
                "${summary_file}" >> ${gene_outfile}
        done
    done

    echo "[DONE] ${g}"

done < ${gene_list}

echo "Master file written to: ${fimo_output_masterfile}"

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



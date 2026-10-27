#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=1000]"
#BSUB -q long
#BSUB -W 08:00
#BSUB -J plot7_submitter
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNB_and_TNF/logs/plot7_submitter_%J.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNB_and_TNF/logs/plot7_submitter_%J.err"
# Plot every eGene lead-variant pair (eigenMT FDR<0.05) as a job array, then
# combine into one PDF. Uses the step5 genotype + CPM SUBSET (which already
# contains ALL donors for the eGene set) so each array task loads a small table
# instead of full chunks. Run step5_get_genotype_and_CPM_of_eGenes.sh first.
set -e
ROOT=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNB_and_TNF
STIMS=IFNB,TNF
SCRIPTDIR=${ROOT}/scripts
RES=${ROOT}/results
GENO=${RES}/genotype_of_eGene_tested_variants.txt   # from step5_get_genotype_and_CPM_of_eGenes.sh
CPM=${RES}/CPM_of_eGenes.txt                          # (both contain all donors/samples)
STATS=${RES}/eGene_variant_stats_all_contexts.txt
PLOTDIR=${ROOT}/plots
mkdir -p ${PLOTDIR}/temp_output ${ROOT}/logs

for f in ${GENO} ${CPM}; do
  [ -s "$f" ] || { echo "ERROR: missing ${f} -- run step5_get_genotype_and_CPM_of_eGenes.sh first." >&2; exit 1; }
done

# ---- combined LONG stats across all contexts (one header) ----
first_stats=$(ls ${RES}/QTL_output_of_all_tested_variants_of_eGenes_*.txt 2>/dev/null | head -1)
[ -n "${first_stats}" ] || { echo "ERROR: no step5 compiled eGene outputs in ${RES}." >&2; exit 1; }
head -1 ${first_stats} > ${STATS}
for f in ${RES}/QTL_output_of_all_tested_variants_of_eGenes_*.txt; do awk 'FNR>1' ${f} >> ${STATS}; done

# ---- unique eGene lead-variant pairs from eigenMT FDR<0.05 tables ----
pairs=${PLOTDIR}/egene_lead_pairs.txt
: > ${pairs}
for f in ${RES}/*_gene_fdr05_table.txt; do
  [ -e "$f" ] || continue
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="gene")gc=i; if($i=="lead_snp")lc=i}; next}
              (gc && lc && $lc!="" && $lc!="NA"){print $gc"\t"$lc}' ${f} >> ${pairs}
done
sort -u ${pairs} -o ${pairs}

# ---- build a self-contained command per pair (one line each) ----
# Writing full commands to a file (and eval-ing one per array index) avoids the
# fragile multi-line `bash -lc` quoting that LSF mangles in the job spool.
CONTAINER=/share/pkg/containers/rstudio_example/r_4.5.2.sif
RSCRIPT=${SCRIPTDIR}/step7_make_plot_gene_snp_pair.R
cmds=${PLOTDIR}/plot7_commands.txt
: > ${cmds}
while IFS=$'\t' read -r g snp; do
  [ -z "${g}" ] && continue
  OUT="${PLOTDIR}/temp_output/plot_${g}_${snp}.pdf"
  echo "singularity exec ${CONTAINER} Rscript ${RSCRIPT} ${ROOT} ${g} ${snp} ${GENO} ${CPM} ${STATS} ${OUT} ${STIMS}" >> ${cmds}
done < ${pairs}
N=$(wc -l < ${cmds})
[ "${N}" -gt 0 ] || { echo "ERROR: no eGene-lead pairs found in ${RES}/*_gene_fdr05_table.txt" >&2; exit 1; }
echo "unique eGene-lead pairs: ${N}"

# ---- submit the plotting array (one command line per index) ----
bsub -q short -W 0:30 -n 1 -R "span[hosts=1]" -R "rusage[mem=3000]" \
  -J "plot7[1-${N}]%200" \
  -o "${ROOT}/logs/plot7_%J_%I.out" -e "${ROOT}/logs/plot7_%J_%I.err" \
  "export R_LIBS_USER=\$HOME/R/x86_64-pc-linux-gnu-library/4.2; cmd=\$(sed -n \"\${LSB_JOBINDEX}p\" ${cmds}); echo \"RUN: \$cmd\"; eval \"\$cmd\""

# ---- wait for the array, then combine into one PDF ----
sleep 15
while [[ $(bjobs -w 2>/dev/null | grep -c "plot7\[") -ne 0 ]]; do
  echo "$(bjobs -w 2>/dev/null | grep -c 'plot7\[') plot jobs remaining"; sleep 15
done
cd ${PLOTDIR}/temp_output
module load poppler 2>/dev/null || true
if command -v pdfunite >/dev/null 2>&1; then
  pdfunite plot_*.pdf ${PLOTDIR}/all_eGene_lead_pairs.pdf && echo "combined -> ${PLOTDIR}/all_eGene_lead_pairs.pdf"
else
  echo "pdfunite not available; individual PDFs left in ${PLOTDIR}/temp_output/"
fi

#!/bin/bash
# submit QTL pipeline steps -- MATCHED-DONOR variant (donors with IFNB & TNF expr + genotype)
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=1000]"
#BSUB -W 72:00
#BSUB -q long
#BSUB -J pipeline
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNB_and_TNF/logs/pipeline_QTL_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNB_and_TNF/logs/pipeline_QTL_%J_%I.err"

#
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/donors_with_IFNB_and_TNF
SCRIPTDIR=${DIR}/scripts
RESULTDIR=${DIR}/results

# ---- PREREQUISITES (run once, in order) ----
# NOTE: this pipeline builds its OWN chunks over the union of IFNB/TNF
# expressed-&-induced genes (25 genes/chunk). It reuses only the genome-wide
# gene-TSS bed and the per-context induced gene lists (read-only), so the
# genome-wide step0.1 must have already produced those.
#
# 0) Build output tree + donor subset + union gene set + gene/pairs/genotype/VST
#    chunks. This step prints "Num gene chunks: N" -- set the [1-N] array upper
#    bound in step1/step2/step3 submit scripts to that N before submitting.
#bsub < ${SCRIPTDIR}/step0_preQTLPipelineProcessing.sh

# 1) Fit QTL models on the matched-donor subset (N-chunk job array):
bsub < ${SCRIPTDIR}/step1_submit_QTL_chunk.sh

#   sanity check for failed chunks (edit job id):
cd ${DIR}/FRB/results   # (repeat per celltype)
for f in ${DIR}/logs/qtl_*.out; do
  if ! grep -q "chunk done" "$f"; then echo "$f"; fi
done

bsub < ${SCRIPTDIR}/step2_QTL_result_QC.sh
echo "done with result QC!"; date

echo "sleeping now for 30min"; date
sleep 1800  # wait for 30min before running eigenMT
bsub < ${SCRIPTDIR}/step3_QTL_result_correction_eigenMT.sh
echo "done with eigenMT!"; date

echo "sleeping now for 30min"; date
sleep 1800 # wait for 30min before running BH correction
bsub < ${SCRIPTDIR}/step4_QTL_result_correction_BH.sh
echo "done with BH!"; date


cd ${RESULTDIR}
cat eGene_summary_by_context.txt
for f in *e_BH_table.txt;do head -2 ${f} | tail -1 | cut -f1-3,7 >> donor_size_summary_by_context.txt; done

# compile eGene and non-eGene and their tested variants that passed QC
bsub < ${SCRIPTDIR}/step5_compile_QTL_output_of_eGenes.sh
bsub < ${SCRIPTDIR}/step6_match_control_nonEgenes.sh

# when last step is done running, also get CPM and genotype of the eGenes and matched-non-eGenes and their QC-passed tested variants
bsub < ${SCRIPTDIR}/step5_get_genotype_and_CPM_of_eGenes.sh
bsub < ${SCRIPTDIR}/step6_get_genotype_and_CPM_of_nonEgenes.sh

# when last step is done, make plots
bsub < ${SCRIPTDIR}/step7_make_plot_all_eGene_lead_pairs.sh


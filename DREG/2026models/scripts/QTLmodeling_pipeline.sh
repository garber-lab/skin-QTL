#!/bin/bash
# submit QTL pipeline steps
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=1000]"
#BSUB -W 72:00
#BSUB -q long
#BSUB -J pipeline
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/pipeline_QTL_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models/logs/pipeline_QTL_%J_%I.err"


# 
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models
SCRIPTDIR=${DIR}/scripts
DATADIR=${DIR}/data
RESULTDIR=${DIR}/results

#bsub < ${SCRIPTDIR}/step0_preQTLPipelineProcessing.sh
#
#bsub < ${SCRIPTDIR}/step1_submit_QTL_chunk.sh #
#
#cd ${DATADIR}/chunks
#awk 'FNR==1 {next} {print $4 "\t" FILENAME}' pairs_chunk_*.tsv \
#  | sed 's/pairs_chunk_//; s/\.tsv//' \
#  | sort -u > gene_chunk_lookup_table.tsv
#
#for f in qtl_428994_*.out; do
#  if ! grep -q "Successfully" "$f"; then
#    echo "$f"
#  fi
#done

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

bsub < ${SCRIPTDIR}/step5_compile_QTL_output_of_eGenes.sh
bsub < ${SCRIPTDIR}/step6_match_control_nonEgenes.sh

bsub < ${SCRIPTDIR}/step5_get_genotype_and_CPM_of_eGenes.sh
bsub < ${SCRIPTDIR}/step6_get_genotype_and_CPM_of_nonEgenes.sh

bsub < ${SCRIPTDIR}/step7_make_plot_all_eGene_lead_pairs.sh

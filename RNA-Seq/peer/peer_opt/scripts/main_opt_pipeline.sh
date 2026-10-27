# steps to set-up and run analysis: optimizing the number of covariates for QTL modeling for each celltype

SIF=/share/pkg/containers/rstudio_example/r_4.5.2.sif
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer
SCRIPTDIR=${DIR}/peer_opt/scripts
PROD=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models
VST=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt

# --- 1. Smoke-test ONE setting end-to-end before the sweep ---
# edit the CONFIG block at the top of opt_submit.sh first, then:
# build a 5-gene, 1-sub-chunk smoke fold (~20 min to fit)
singularity exec $SIF Rscript ${SCRIPTDIR}/opt_make_folds.R \
  --chunks-dir $PROD/data/chunks --out-dir ${DIR}/peer_opt/folds_smoke/FRB \
  --n-folds 1 --n-genes-per-fold 5 \
  --seed-genes-file $PROD/results/FRB_PBS_eQTL_gene_fdr05_genelist.txt \
  --n-seed-per-fold 2 --subchunk-genes 5 --vst $VST --seed 42

unset LSB_JOBINDEX; SMOKE=1 bash ${SCRIPTDIR}/opt_submit_smoke.sh
# runs combo #1 locally (no bsub) and appends one row-group to the summary.


# --- 2. build folds --- #
# per celltype, 500 genes, 3 folds. 120/500 gene are known eGenes. 380 genes are random genes selected from the rest of the gene list. Each fold is sub-chunked into 25 gene chunks.
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
for ct in FRB KRT MEL; do
  singularity exec $SIF Rscript ${SCRIPTDIR}/opt_make_folds.R \
    --chunks-dir $PROD/data/chunks --out-dir ${DIR}/peer_opt/folds/$ct \
    --n-folds 3 --n-genes-per-fold 500 \
    --seed-genes-file $PROD/results/${ct}_PBS_eQTL_gene_fdr05_genelist.txt \
    --n-seed-per-fold 120 --subchunk-genes 25 --vst $VST --seed 42
done

# --- 3. Submit the full grid --- #
unset LSB_JOBINDEX; bash ${SCRIPTDIR}/opt_submit.sh          # bsub array [1-540]%192

# 4. Summary



DIR=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

# --- generate a grid of expression PEER or expression PCs --- #
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
Rscript ${DIR}/generate_covariates.R \
  --expression /pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt \
  --metadata   /pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt \
  --cpm        /pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt \
  --min-cpm 1 --min-samples 5 \
  --all-cell-types --k-grid 0,5,8,10,15 --method both \
  --output-dir ${DIR}/peer_opt --seed 42

# --- also create expression PCs with stimulation condition regressed out to prevent colinearity --- #
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
Rscript ${DIR}/generate_covariates.R \
  --expression /pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt \
  --metadata   /pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt \
  --cpm        /pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt \
  --min-cpm 1 --min-samples 5 \
  --all-cell-types --k-grid 0,5,8,10,15 --method expression_pc --pc-remove-condition \
  --output-dir ${DIR}/peer_opt --seed 42

# --- annotate each PEER or expression PC and create figures --- #
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
Rscript ${DIR}/compare_covariate_annotation.R \
  --output-dir ${DIR}/peer_opt \
  --metadata /pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt \
  --vars condition,sequencing_platform,phase,read_length,total_reads,unique_reads_aligned,all_aligned_reads_over_total,CS2 \
  2>&1 | grep -Ev "Attaching|masked|following|%notin%|^$|package:" | sed 's/^/  /'
echo; echo "=== outputs ==="
ls -1 ${DIR}/peer_opt/covariates/annotation_compare_*.pdf

#=== outputs ===
#/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_opt/covariates/annotation_compare_FRB_expression_pc_condremoved.pdf
#/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_opt/covariates/annotation_compare_FRB_expression_pc.pdf
#/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_opt/covariates/annotation_compare_FRB_peer.pdf
#/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_opt/covariates/annotation_compare_KRT_expression_pc_condremoved.pdf
#/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_opt/covariates/annotation_compare_KRT_expression_pc.pdf
#/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_opt/covariates/annotation_compare_KRT_peer.pdf
#/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_opt/covariates/annotation_compare_MEL_expression_pc_condremoved.pdf
#/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_opt/covariates/annotation_compare_MEL_expression_pc.pdf
#/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_opt/covariates/annotation_compare_MEL_peer.pdf


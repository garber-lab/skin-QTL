# pipeline to re-run the QTL modeling pipeline but only for donors that are present in both IFNB and IFNG

celltype=FRB
workingdir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${celltype}/IFNBeGene_depletion_issue
mkdir -p ${workingdir}
mkdir -p ${workingdir}/data ${workingdir}/logs ${workingdir}/chunks ${workingdir}/results

# --- 1. only restricting to induced genes in either IFNB and IFNG --- #
#/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_induced_genes_log2FC1_padj0.05.R
#fetch_induced_genes(this_celltype, "IFNG")
# and fetch_expressed_genes(this_celltype, "IFNG")
#already written to /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/IFNBeGene_depletion_issue/data/induced_genes_FRB_IFNB.txt 
# ct="MEL"; condition="TNF"; expressed=fetch_expressed_genes(ct, condition); induced=fetch_induced_genes(ct, condition); induced_and_expressed=intersect(induced, expressed); write(induced_and_expressed, paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct,"/IFNBeGene_depletion_issue/data/expressed_and_induced_genes_",ct,"_",condition,".txt"))
cat ${workingdir}/data/expressed_and_induced_genes_${celltype}_IFNB.txt ${workingdir}/data/expressed_and_induced_genes_${celltype}_IFNG.txt | sort -u > ${workingdir}/data/expressed_and_induced_genes_${celltype}_IFNB_or_IFNG.txt


# --- 2. chunk gene-variant pairs --- #
bash step16_IFNBQTLdepletion_QTLmodel_preQTLPipelineProcessing.sh
# there are 119 chunks with max 20,000 pairs each chunk

# --- 3. submit pipeline across chunks --- #
# takes about 4 hours to run
bsub < step16_IFNBQTLdepletion_QTLmodel_submit_chunks.sh
mkdir -p ${workingdir}/logs/qtl; mv ${workingdir}/logs/qtl_* ${workingdir}/logs/qtl


# --- 4. post QTLmodeling QC --- #
bsub < step16_IFNBQTLdepletion_QTLmodel_QTLresultQC.sh
mkdir -p ${workingdir}/logs/qc; mv ${workingdir}/logs/qc_* mv ${workingdir}/logs/qc

# --- 5. multiple testing correction --- #
bsub < step16_IFNBQTLdepletion_QTLmodel_eigenMT_perGeneChunk.sh
bsub < step16_IFNBQTLdepletion_QTLmodel_BHcorrection.sh 
mkdir -p ${workingdir}/logs/eigenMT; mv ${workingdir}/logs/eigenMT_* mv ${workingdir}/logs/eigenMT
mkdir -p ${workingdir}/logs/BH; mv ${workingdir}/logs/BH_* mv ${workingdir}/logs/BH


# --- 6. analyze IFNB depletion via statistical test --- #
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
  Rscript step16_IFNBQTLdepletion_QTLmodel_McNemarTest.R    # no args = FRB IFNG IFNB eQTL 0.05

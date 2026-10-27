#!/bin/bash
module load conda; source activate fastQTL
set -euo pipefail

# motif scan (sequence only, model-independent); a FIMO failure should not block the tables/plots
bash /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/step6_fimo_snp.sh /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/posterchild/rs9676826 50 || echo "WARNING: FIMO step failed, continuing" >&2
singularity exec --nv --env PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages --env TF_FORCE_GPU_ALLOW_GROWTH=true /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif bash -c "export PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages; python /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/step6_summarize_snp_contrib.py     --outdir /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/posterchild/rs9676826 --rsid rs9676826 --cb /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet --bias 0.5     --models FRB_IFN,FRB_PBS,KRT_IFN,KRT_PBS,MEL_IFN,MEL_PBS --folds 0,1,2,3,4 --genome /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/data/downloads/hg38.fa"

# plots: one PDF per head -> /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/posterchild/rs9676826/rs9676826.contrib_<head>.pdf
export R_LIBS_USER=/home/shuo.shan-umw/R/x86_64-pc-linux-gnu-library/4.2

# motif annotation: Fi-NeMo hits + MoDISco names + FIMO + contributions/ISM + expression
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif Rscript /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/step6_motif_annotate.R /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/posterchild/rs9676826 /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet 0.5 200 50   || echo "WARNING: motif annotation failed, plotting without the motif row" >&2

for h in counts profile; do
  singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif Rscript /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/step6_plot_snp_contrib.R /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/posterchild/rs9676826 ${h} 45 250
done

# one merged PDF: ChromBPNet figures first, then the FIMO report per context
if [[ "1" == "1" ]]; then
  pdfs=()
  for h in counts profile; do
    f=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/posterchild/rs9676826/rs9676826.contrib_${h}.pdf; if [[ -s ${f} ]]; then pdfs+=("${f}"); fi
  done
  for ctx in KRT_PBS KRT_IFNG FRB_PBS FRB_IFNG MEL_PBS MEL_IFNG; do
    f=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/posterchild/rs9676826/motifs/fimo/fimo_summary_rs9676826_${ctx}.pdf; if [[ -s ${f} ]]; then pdfs+=("${f}"); fi
  done
  merged=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/posterchild/rs9676826/rs9676826.report.pdf
  if [[ ${#pdfs[@]} -eq 0 ]]; then
    echo "WARNING: no PDFs to merge" >&2
  else
    # R qpdf: the only merger available on this cluster (no pdfunite/gs/qpdf binaries)
    singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif Rscript -e 'a <- commandArgs(TRUE); qpdf::pdf_combine(a[-1], a[1])'       "${merged}" "${pdfs[@]}" && echo "merged ${#pdfs[@]} PDFs -> ${merged}"
  fi
fi

#!/bin/bash
module load conda; source activate fastQTL
set -euo pipefail
rsid=$(cut -f5 /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/snpList/IFNG_shortlist_variants_rsid/variants.tsv | sed -n "${LSB_JOBINDEX}p")
[[ -n "${rsid}" ]] || { echo "ERROR: no rsID at index ${LSB_JOBINDEX}" >&2; exit 1; }
vdir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/posterchild/${rsid}
echo "== ${rsid} -> ${vdir}"

bash /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/step6_fimo_snp.sh ${vdir} 50 || echo "WARNING: FIMO step failed for ${rsid}" >&2

singularity exec --nv --env PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages --env TF_FORCE_GPU_ALLOW_GROWTH=true /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif bash -c "export PYTHONPATH=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages; python /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/step6_summarize_snp_contrib.py     --outdir ${vdir} --rsid ${rsid} --cb /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet --bias 0.5     --models FRB_IFN,FRB_PBS,KRT_IFN,KRT_PBS,MEL_IFN,MEL_PBS --folds 0,1,2,3,4 --genome /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/data/downloads/hg38.fa"

export R_LIBS_USER=/home/shuo.shan-umw/R/x86_64-pc-linux-gnu-library/4.2
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif Rscript /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/step6_motif_annotate.R ${vdir} /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet 0.5 200 50   || echo "WARNING: motif annotation failed for ${rsid}" >&2
for h in counts profile; do
  singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif Rscript /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/step6_plot_snp_contrib.R ${vdir} ${h} 45 250
done

if [[ "1" == "1" ]]; then
  pdfs=()
  for h in counts profile; do
    f=${vdir}/${rsid}.contrib_${h}.pdf; if [[ -s ${f} ]]; then pdfs+=("${f}"); fi
  done
  for ctx in KRT_PBS KRT_IFNG FRB_PBS FRB_IFNG MEL_PBS MEL_IFNG; do
    f=${vdir}/motifs/fimo/fimo_summary_${rsid}_${ctx}.pdf; if [[ -s ${f} ]]; then pdfs+=("${f}"); fi
  done
  if [[ ${#pdfs[@]} -gt 0 ]]; then
    singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif Rscript -e 'a <- commandArgs(TRUE); qpdf::pdf_combine(a[-1], a[1])'       "${vdir}/${rsid}.report.pdf" "${pdfs[@]}" && echo "merged ${#pdfs[@]} PDFs"
  fi
fi
echo "done ${rsid}"

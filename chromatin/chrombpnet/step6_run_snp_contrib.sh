#!/bin/bash
# Per-variant ChromBPNet contribution scores + allele predictions
# Usage:  bash step6_run_snp_contrib.sh <rsID> [outdir]
# Dry run (write job scripts, submit nothing):  DRYRUN=1 bash step6_run_snp_contrib.sh <rsID>
# Output: ${CB}/posterchild/<rsID>/ (all 6 models x 5 folds for that variant)
set -euo pipefail

rsid=$1
DRYRUN=${DRYRUN:-0}
CB=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet
outdir=${2:-${CB}/posterchild/${rsid}}   # one folder per variant; all models inside

############################ CONFIG ############################
MODELS="FRB_IFN FRB_PBS KRT_IFN KRT_PBS MEL_IFN MEL_PBS"
FOLDS="0 1 2 3 4"
BIAS=0.5
VCF=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz
DL=${CB}/ATAC_KRT_IFN_merged_alldonors/biasFactor${BIAS}/data/downloads
GENOME=${DL}/hg38.fa
CHROMSIZES=${DL}/hg38.chrom.sizes
SIF=${CB}/chrombpnet_latest.sif
VS=${CB}/variant-scorer/src
SUMMARIZE_PY=${CB}/step6_summarize_snp_contrib.py   # full path to the summarize script
ISM_PY=${CB}/step6_ism_snp.py                       # full path to the ISM script
ISM_FLANK=50                                         # ISM window: +/- bp around the variant
FIMO_SH=${CB}/step6_fimo_snp.sh                      # FIMO REF vs ALT (HOCOMOCO v12), host-side
FIMO_FLANK=50                                        # FIMO window: +/- bp around the variant
FINEMO=/home/shuo.shan-umw/.conda/envs/finemo/bin/finemo   # Fi-NeMo (own conda env, PyTorch)
# Benchmarked on rs9604525 (KRT_IFN, profile): w=200 loose gives the same hits over +/-50 bp as
# w=400 tight in ~100 s instead of ~233 s. w=100 is too narrow: it truncates CWMs and scrambles the fit.
FINEMO_WIDTH=200                                     # hit-calling window (variant at its centre)
FINEMO_OPTS="-c 0.0001 -S 10000"                     # convergence tolerance / max steps
FINEMO_DEVICE=cpu                                    # cpu sets CUDA_VISIBLE_DEVICES="": reliable when
                                                     # two model jobs share a GPU with TensorFlow, and
                                                     # deterministic. gpu is faster but never mix the
                                                     # two devices within one analysis
FINEMO_KEEP_FLANK=50                                 # annotate only hits within +/- this of the variant
FINEMO_HEADS="counts profile"                        # contribution heads to call hits on
# motif library per model: pipeline MoDISco (fold0, profile head, 30K peaks) -> swap for the
# full-peak fold-averaged MoDISco later; glob is relative to each model's biasFactor folder
MODISCO_GLOB="fold0/chrombpnet_model/auxiliary/interpret_subsample/*modisco_results_profile_scores.h5"

# Plotting (runs at the end of the summary job, in the R container)
PLOT_R=${CB}/step6_plot_snp_contrib.R                      # full path to the plotting script
MOTIF_R=${CB}/step6_motif_annotate.R                 # joins Fi-NeMo + FIMO + contributions + expression
R_SIF=/share/pkg/containers/rstudio_example/r_4.5.2.sif
R_LIBS=${HOME}/R/x86_64-pc-linux-gnu-library/4.2    # needs ggseqlogo, patchwork, gridExtra, tidyverse
PLOT_HEADS="counts profile"                          # one PDF per contribution head
LOGO_FLANK=45                                        # logo / ISM heatmap / motif window (<= ISM_FLANK)
TRACK_FLANK=250                                      # observed / predicted track window
MERGE_PDF=1                                          # 1 = also write <rsID>.report.pdf (all figures in one)
MERGE_ORDER="KRT_PBS KRT_IFNG FRB_PBS FRB_IFNG MEL_PBS MEL_IFNG"   # FIMO report order after the figures

# Environment copied from step3_variant_scoring_KRT_IFNG_merged.sh (known working)
CONDA_SITE=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages
HOST_SETUP="module load conda; source activate fastQTL"        # run on the node before singularity
# TF_FORCE_GPU_ALLOW_GROWTH: TensorFlow otherwise pre-allocates ~31 GB of the 32 GB V100, so a second
# model job landing on the same node crashes with CUDA_ERROR_ILLEGAL_ADDRESS (exit 134)
RUN="singularity exec --nv --env PYTHONPATH=${CONDA_SITE} --env TF_FORCE_GPU_ALLOW_GROWTH=true ${SIF}"
ENV_SETUP="export PYTHONPATH=${CONDA_SITE};"                     # run inside the container
# Sized from the rs9604525 runs: max 2.8 GB RAM; ~13 min per model job with ISM + Fi-NeMo on CPU,
# less on GPU. mem is per core here: 2 cores x 3000 MB = 6 GB (~2x peak).
BSUB_OPTS='-q gpu -W 1:00 -n 2 -R "rusage[mem=3000]" -R "span[hosts=1]"'
################################################################

if [[ ! -f "${SUMMARIZE_PY}" ]]; then
  echo "ERROR: summarize script not found: ${SUMMARIZE_PY} (set SUMMARIZE_PY in the config)" >&2; exit 1
fi
if [[ ! -f "${MOTIF_R}" ]]; then
  echo "ERROR: motif annotation script not found: ${MOTIF_R} (set MOTIF_R in the config)" >&2; exit 1
fi
if [[ ! -f "${PLOT_R}" ]]; then
  echo "ERROR: plot script not found: ${PLOT_R} (set PLOT_R in the config)" >&2; exit 1
fi
if [[ ! -x "${FINEMO}" ]]; then
  echo "ERROR: finemo not found/executable: ${FINEMO} (set FINEMO in the config)" >&2; exit 1
fi
if [[ ! -f "${FIMO_SH}" ]]; then
  echo "ERROR: FIMO script not found: ${FIMO_SH} (set FIMO_SH in the config)" >&2; exit 1
fi
if [[ ! -f "${ISM_PY}" ]]; then
  echo "ERROR: ISM script not found: ${ISM_PY} (set ISM_PY in the config)" >&2; exit 1
fi

# device switch for Fi-NeMo (empty CUDA_VISIBLE_DEVICES hides the GPU)
if [[ "${FINEMO_DEVICE}" == "cpu" ]]; then FINEMO_ENV='CUDA_VISIBLE_DEVICES=""'; else FINEMO_ENV=""; fi

module load bcftools

mkdir -p "${outdir}/per_fold" "${outdir}/logs" "${outdir}/jobs"
cd "${outdir}"

### 1. variant file, chrombpnet schema, no header:
###    chr  pos(1-based)  allele1(REF)  allele2(ALT)  variant_id
varF=${outdir}/${rsid}.variant.tsv
bcftools view -i "ID=\"${rsid}\"" "${VCF}" \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' > "${varF}"

n=$(wc -l < "${varF}")
if [[ ${n} -ne 1 ]]; then
  echo "ERROR: expected 1 VCF record for ${rsid}, found ${n}" >&2; cat "${varF}" >&2; exit 1
fi
if cut -f4 "${varF}" | grep -q ','; then
  echo "ERROR: ${rsid} is multi-allelic; split it first" >&2; exit 1
fi
echo "Variant:"; cat "${varF}"

### 2. one job per model; folds run sequentially inside the job
job_names=()
for m in ${MODELS}; do
  base=${CB}/ATAC_${m}_merged_alldonors/biasFactor${BIAS}
  jobF=${outdir}/jobs/${m}.sh
  cat > "${jobF}" <<EOF
#!/bin/bash
${HOST_SETUP}
set -euo pipefail
shopt -s nullglob
cd ${VS}
for f in ${FOLDS}; do
  # pipeline prefixes model files with <sample>_fold<f>_, so glob for the nobias model
  hits=(${base}/fold\${f}/chrombpnet_model/models/*chrombpnet_nobias.h5)
  if [[ \${#hits[@]} -ne 1 ]]; then
    echo "ERROR: expected 1 *chrombpnet_nobias.h5 in fold\${f}, found \${#hits[@]}: \${hits[*]}" >&2
    exit 1
  fi
  model=\${hits[0]}
  pre=${outdir}/per_fold/${m}_fold\${f}
  echo "fold\${f} model: \${model}"

  # REF/ALT contribution scores (counts + profile heads)
  ${RUN} bash -c "${ENV_SETUP} python variant_shap.py \
      -l ${varF} -g ${GENOME} -m \${model} -o \${pre} -s ${CHROMSIZES} \
      -sc chrombpnet -st counts profile -bs 1"

  # REF/ALT predicted profiles + counts (saved to \${pre}.variant_predictions.h5)
  ${RUN} bash -c "${ENV_SETUP} python variant_scoring.py \
      -l ${varF} -g ${GENOME} -m \${model} -o \${pre} -s ${CHROMSIZES} \
      -sc chrombpnet -n 10"
done

# saturation mutagenesis around the variant, all folds in one python process
${RUN} bash -c "${ENV_SETUP} python ${ISM_PY} --vs_src ${VS} \
    --variant ${varF} --genome ${GENOME} \
    --model_glob '${base}/fold{f}/chrombpnet_model/models/*chrombpnet_nobias.h5' \
    --folds ${FOLDS// /,} --flank ${ISM_FLANK} --out ${outdir}/per_fold/${m}.ism.tsv"

# Fi-NeMo motif hits on fold-averaged REF/ALT contributions (region 0 = REF, 1 = ALT),
# against this model's own MoDISco patterns; failures here should not fail the job
mh=(${base}/${MODISCO_GLOB})
if [[ \${#mh[@]} -ne 1 ]]; then
  echo "WARNING: expected 1 MoDISco h5 for ${m}, found \${#mh[@]}; skipping Fi-NeMo" >&2
else
  for h in ${FINEMO_HEADS}; do
    fd=${outdir}/motifs/finemo/${m}/\${h}
    mkdir -p \${fd}
    h5s=(${outdir}/per_fold/${m}_fold*.variant_shap.\${h}.h5)
    { ${FINEMO_ENV} ${FINEMO} extract-regions-chrombpnet-h5 -c "\${h5s[@]}" -w ${FINEMO_WIDTH} -o \${fd}/regions.npz &&
      ${FINEMO_ENV} ${FINEMO} call-hits -r \${fd}/regions.npz -m \${mh[0]} -o \${fd} ${FINEMO_OPTS} &&
      [[ -s \${fd}/hits.tsv ]] &&
      echo "\${mh[0]}" > \${fd}/modisco_source.txt; } \
      || echo "WARNING: Fi-NeMo failed for ${m} \${h}" >&2
  done
fi
echo "done ${m}"
EOF
  jname=snpc_${rsid}_${m}
  if [[ ${DRYRUN} -eq 1 ]]; then echo "[dry run] wrote ${jobF}"; continue; fi
  # capture the job ID: name-based dependencies break once LSF purges finished jobs
  sub=$(eval bsub ${BSUB_OPTS} -J "${jname}" -oo "${outdir}/logs/${m}.out" -eo "${outdir}/logs/${m}.err" bash "${jobF}")
  echo "${sub}"
  jid=$(sed -n 's/^Job <\([0-9]*\)>.*/\1/p' <<< "${sub}")
  [[ -n "${jid}" ]] || { echo "ERROR: could not parse job id from: ${sub}" >&2; exit 1; }
  job_names+=("ended(${jid})")   # ended, not done: a failed model job should not strand the summary
done

### 3. summary job, runs after all model jobs finish
sumF=${outdir}/jobs/summary.sh
cat > "${sumF}" <<EOF
#!/bin/bash
${HOST_SETUP}
set -euo pipefail

# motif scan (sequence only, model-independent); a FIMO failure should not block the tables/plots
bash ${FIMO_SH} ${outdir} ${FIMO_FLANK} || echo "WARNING: FIMO step failed, continuing" >&2
${RUN} bash -c "${ENV_SETUP} python ${SUMMARIZE_PY} \
    --outdir ${outdir} --rsid ${rsid} --cb ${CB} --bias ${BIAS} \
    --models ${MODELS// /,} --folds ${FOLDS// /,} --genome ${GENOME}"

# plots: one PDF per head -> ${outdir}/${rsid}.contrib_<head>.pdf
export R_LIBS_USER=${R_LIBS}

# motif annotation: Fi-NeMo hits + MoDISco names + FIMO + contributions/ISM + expression
singularity exec ${R_SIF} Rscript ${MOTIF_R} ${outdir} ${CB} ${BIAS} ${FINEMO_WIDTH} ${FINEMO_KEEP_FLANK} \
  || echo "WARNING: motif annotation failed, plotting without the motif row" >&2

for h in ${PLOT_HEADS}; do
  singularity exec ${R_SIF} Rscript ${PLOT_R} ${outdir} \${h} ${LOGO_FLANK} ${TRACK_FLANK}
done

# one merged PDF: ChromBPNet figures first, then the FIMO report per context
if [[ "${MERGE_PDF}" == "1" ]]; then
  pdfs=()
  for h in ${PLOT_HEADS}; do
    f=${outdir}/${rsid}.contrib_\${h}.pdf; if [[ -s \${f} ]]; then pdfs+=("\${f}"); fi
  done
  for ctx in ${MERGE_ORDER}; do
    f=${outdir}/motifs/fimo/fimo_summary_${rsid}_\${ctx}.pdf; if [[ -s \${f} ]]; then pdfs+=("\${f}"); fi
  done
  merged=${outdir}/${rsid}.report.pdf
  if [[ \${#pdfs[@]} -eq 0 ]]; then
    echo "WARNING: no PDFs to merge" >&2
  else
    # R qpdf: the only merger available on this cluster (no pdfunite/gs/qpdf binaries)
    singularity exec ${R_SIF} Rscript -e 'a <- commandArgs(TRUE); qpdf::pdf_combine(a[-1], a[1])' \
      "\${merged}" "\${pdfs[@]}" && echo "merged \${#pdfs[@]} PDFs -> \${merged}"
  fi
fi
EOF
if [[ ${DRYRUN} -eq 1 ]]; then echo "[dry run] wrote ${sumF}"; exit 0; fi

dep=$(printf ' && %s' "${job_names[@]}"); dep=${dep:4}
bsub -q gpu -W 0:30 -n 1 -R "rusage[mem=4000]" \
  -J "snpc_${rsid}_summary" -w "${dep}" \
  -oo "${outdir}/logs/summary.out" -eo "${outdir}/logs/summary.err" \
  bash "${sumF}"

echo "Submitted. Results -> ${outdir}"

#!/bin/bash
# ChromBPNet per-variant analysis for a LIST of variants (batched version of step6_run_snp_contrib.sh)
#
# Usage:  bash step6b_snpList_run.sh <rsid_list.txt> [batch_name]
# Dry run (write job scripts, submit nothing):  DRYRUN=1 bash step6b_snpList_run.sh <rsid_list.txt>
#
# Stage 1  build one variant table for the whole list (single bcftools pass)
# Stage 2  6 jobs, one per model: variant_shap + variant_scoring (all variants, per fold),
#          batched ISM (one model load per fold), Fi-NeMo on all 2N regions at once
# Stage 3  split the batched outputs into ${CB}/posterchild/<rsID>/ (step6 layout)
# Stage 4  LSF array, one index per variant: FIMO (+ reports), summarize, motif annotation,
#          both figures, merged PDF - all by calling the unchanged step6 scripts
set -euo pipefail

LIST=$(readlink -f "$1")
BATCH_NAME=${2:-$(basename "${LIST}" .txt)}
DRYRUN=${DRYRUN:-0}
CB=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet

############################ CONFIG ############################
MODELS="FRB_IFN FRB_PBS KRT_IFN KRT_PBS MEL_IFN MEL_PBS"
FOLDS="0 1 2 3 4"
BIAS=0.5
OUTROOT=${CB}/posterchild                            # per-variant results (step6 layout)
BATCH=${CB}/snpList/${BATCH_NAME}                    # batched intermediates + logs
VCF=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz
DL=${CB}/ATAC_KRT_IFN_merged_alldonors/biasFactor${BIAS}/data/downloads
GENOME=${DL}/hg38.fa
CHROMSIZES=${DL}/hg38.chrom.sizes
SIF=${CB}/chrombpnet_latest.sif
VS=${CB}/variant-scorer/src

# step6b scripts (batched) and step6 scripts (per variant, unchanged)
ISM_PY=${CB}/step6b_snpList_ism.py
SPLIT_PY=${CB}/step6b_snpList_split.py
SUMMARIZE_PY=${CB}/step6_summarize_snp_contrib.py
FIMO_SH=${CB}/step6_fimo_snp.sh
MOTIF_R=${CB}/step6_motif_annotate.R
PLOT_R=${CB}/step6_plot_snp_contrib.R

SHAP_BS=8                                            # variant_shap batch size (memory/speed)
SCORE_SHUFFLES=10                                    # variant_scoring -n (null shuffles per variant)
ISM_FLANK=50
FIMO_FLANK=50
FINEMO=/home/shuo.shan-umw/.conda/envs/finemo/bin/finemo
FINEMO_WIDTH=200
FINEMO_OPTS="-c 0.0001 -S 10000"
FINEMO_DEVICE=cpu                                    # cpu = deterministic; see step6 notes
FINEMO_HEADS="counts profile"
MODISCO_GLOB="fold0/chrombpnet_model/auxiliary/interpret_subsample/*modisco_results_profile_scores.h5"

R_SIF=/share/pkg/containers/rstudio_example/r_4.5.2.sif
R_LIBS=${HOME}/R/x86_64-pc-linux-gnu-library/4.2
PLOT_HEADS="counts profile"
LOGO_FLANK=45
TRACK_FLANK=250
MERGE_PDF=1
MERGE_ORDER="KRT_PBS KRT_IFNG FRB_PBS FRB_IFNG MEL_PBS MEL_IFNG"

CONDA_SITE=/home/shuo.shan-umw/miniconda3/envs/fastQTL/lib/python3.9/site-packages
HOST_SETUP="module load conda; source activate fastQTL"
RUN="singularity exec --nv --env PYTHONPATH=${CONDA_SITE} --env TF_FORCE_GPU_ALLOW_GROWTH=true ${SIF}"
ENV_SETUP="export PYTHONPATH=${CONDA_SITE};"
# model jobs carry all variants: more time and memory than the single-variant step6 jobs
MODEL_BSUB='-q gpu -W 8:00 -n 4 -R "rusage[mem=8000]" -R "span[hosts=1]"'
SPLIT_BSUB='-q short -W 2:00 -n 1 -R "rusage[mem=16000]"'
ARRAY_BSUB='-q short -W 1:00 -n 2 -R "rusage[mem=4000]" -R "span[hosts=1]"'
ARRAY_CONCURRENCY=20                                 # max simultaneous per-variant jobs
################################################################

for f in "${ISM_PY}" "${SPLIT_PY}" "${SUMMARIZE_PY}" "${FIMO_SH}" "${MOTIF_R}" "${PLOT_R}"; do
  [[ -f "${f}" ]] || { echo "ERROR: missing script ${f}" >&2; exit 1; }
done
[[ -x "${FINEMO}" ]] || { echo "ERROR: finemo not executable: ${FINEMO}" >&2; exit 1; }
[[ -s "${LIST}" ]] || { echo "ERROR: empty rsID list: ${LIST}" >&2; exit 1; }
if [[ "${FINEMO_DEVICE}" == "cpu" ]]; then FINEMO_ENV='CUDA_VISIBLE_DEVICES=""'; else FINEMO_ENV=""; fi

module load bcftools
mkdir -p "${BATCH}/per_fold" "${BATCH}/logs" "${BATCH}/jobs" "${OUTROOT}"

### 1. one variant table for the whole list (single VCF pass)
varF=${BATCH}/variants.tsv
grep -vE '^\s*$' "${LIST}" | tr -d '\r' | sort -u > "${BATCH}/rsids.txt"
n_req=$(wc -l < "${BATCH}/rsids.txt")
bcftools view --include "ID==@${BATCH}/rsids.txt" "${VCF}" \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' > "${varF}"

awk -F'\t' '$4 ~ /,/ {print $5}' "${varF}" > "${BATCH}/multiallelic.txt" || true
if [[ -s "${BATCH}/multiallelic.txt" ]]; then
  echo "ERROR: multi-allelic variants (split them first):" >&2; cat "${BATCH}/multiallelic.txt" >&2; exit 1
fi
cut -f5 "${varF}" | sort -u > "${BATCH}/found.txt"
comm -23 "${BATCH}/rsids.txt" "${BATCH}/found.txt" > "${BATCH}/not_in_vcf.txt" || true
n_found=$(wc -l < "${varF}")
echo "requested ${n_req} rsIDs; found ${n_found} records in the VCF"
[[ -s "${BATCH}/not_in_vcf.txt" ]] && { echo "not found (skipped):"; cat "${BATCH}/not_in_vcf.txt"; }
[[ ${n_found} -gt 0 ]] || { echo "ERROR: no variants found in the VCF" >&2; exit 1; }
awk -F'\t' '{ if (length($3)!=1 || length($4)!=1) print "  non-SNV: "$5" "$3">"$4 }' "${varF}"

### 2. one job per model: all variants, all folds, then batched ISM, then Fi-NeMo on 2N regions
job_deps=()
for m in ${MODELS}; do
  base=${CB}/ATAC_${m}_merged_alldonors/biasFactor${BIAS}
  jobF=${BATCH}/jobs/model_${m}.sh
  cat > "${jobF}" <<EOF
#!/bin/bash
${HOST_SETUP}
set -euo pipefail
shopt -s nullglob
cd ${VS}
for f in ${FOLDS}; do
  hits=(${base}/fold\${f}/chrombpnet_model/models/*chrombpnet_nobias.h5)
  if [[ \${#hits[@]} -ne 1 ]]; then
    echo "ERROR: expected 1 *chrombpnet_nobias.h5 in fold\${f}, found \${#hits[@]}" >&2; exit 1
  fi
  model=\${hits[0]}
  pre=${BATCH}/per_fold/${m}_fold\${f}
  echo "== fold\${f} \${model}"

  ${RUN} bash -c "${ENV_SETUP} python variant_shap.py \
      -l ${varF} -g ${GENOME} -m \${model} -o \${pre} -s ${CHROMSIZES} \
      -sc chrombpnet -st counts profile -bs ${SHAP_BS}"

  ${RUN} bash -c "${ENV_SETUP} python variant_scoring.py \
      -l ${varF} -g ${GENOME} -m \${model} -o \${pre} -s ${CHROMSIZES} \
      -sc chrombpnet -n ${SCORE_SHUFFLES}"
done

# batched ISM: one model load per fold, every variant scored against it
${RUN} bash -c "${ENV_SETUP} python ${ISM_PY} --vs_src ${VS} \
    --variants ${varF} --genome ${GENOME} \
    --model_glob '${base}/fold{f}/chrombpnet_model/models/*chrombpnet_nobias.h5' \
    --folds ${FOLDS// /,} --flank ${ISM_FLANK} --out ${BATCH}/per_fold/${m}.ism.tsv"

# Fi-NeMo: all REF+ALT regions of every variant in one fit per head
mh=(${base}/${MODISCO_GLOB})
if [[ \${#mh[@]} -ne 1 ]]; then
  echo "WARNING: expected 1 MoDISco h5 for ${m}, found \${#mh[@]}; skipping Fi-NeMo" >&2
else
  for h in ${FINEMO_HEADS}; do
    fd=${BATCH}/motifs/finemo/${m}/\${h}
    mkdir -p \${fd}
    h5s=(${BATCH}/per_fold/${m}_fold*.variant_shap.\${h}.h5)
    { ${FINEMO_ENV} ${FINEMO} extract-regions-chrombpnet-h5 -c "\${h5s[@]}" -w ${FINEMO_WIDTH} -o \${fd}/regions.npz &&
      ${FINEMO_ENV} ${FINEMO} call-hits -r \${fd}/regions.npz -m \${mh[0]} -o \${fd} ${FINEMO_OPTS} &&
      [[ -s \${fd}/hits.tsv ]] &&
      echo "\${mh[0]}" > \${fd}/modisco_source.txt; } \
      || echo "WARNING: Fi-NeMo failed for ${m} \${h}" >&2
  done
fi
echo "done ${m}"
EOF
  if [[ ${DRYRUN} -eq 1 ]]; then echo "[dry run] wrote ${jobF}"; continue; fi
  sub=$(eval bsub ${MODEL_BSUB} -J "snpL_${BATCH_NAME}_${m}" \
        -oo "${BATCH}/logs/model_${m}.out" -eo "${BATCH}/logs/model_${m}.err" bash "${jobF}")
  echo "${sub}"
  jid=$(sed -n 's/^Job <\([0-9]*\)>.*/\1/p' <<< "${sub}")
  [[ -n "${jid}" ]] || { echo "ERROR: could not parse job id: ${sub}" >&2; exit 1; }
  job_deps+=("ended(${jid})")
done

### 3. split batched outputs into per-variant folders
splitF=${BATCH}/jobs/split.sh
cat > "${splitF}" <<EOF
#!/bin/bash
${HOST_SETUP}
set -euo pipefail
${RUN} bash -c "${ENV_SETUP} python ${SPLIT_PY} \
    --batch ${BATCH} --outroot ${OUTROOT} --variants ${varF} \
    --models ${MODELS// /,} --folds ${FOLDS// /,} --heads ${FINEMO_HEADS// /,}"
EOF

### 4. per-variant array job: FIMO, summary, motif annotation, figures, merged PDF
arrayF=${BATCH}/jobs/per_variant.sh
cat > "${arrayF}" <<EOF
#!/bin/bash
${HOST_SETUP}
set -euo pipefail
rsid=\$(cut -f5 ${varF} | sed -n "\${LSB_JOBINDEX}p")
[[ -n "\${rsid}" ]] || { echo "ERROR: no rsID at index \${LSB_JOBINDEX}" >&2; exit 1; }
vdir=${OUTROOT}/\${rsid}
echo "== \${rsid} -> \${vdir}"

bash ${FIMO_SH} \${vdir} ${FIMO_FLANK} || echo "WARNING: FIMO step failed for \${rsid}" >&2

${RUN} bash -c "${ENV_SETUP} python ${SUMMARIZE_PY} \
    --outdir \${vdir} --rsid \${rsid} --cb ${CB} --bias ${BIAS} \
    --models ${MODELS// /,} --folds ${FOLDS// /,} --genome ${GENOME}"

export R_LIBS_USER=${R_LIBS}
singularity exec ${R_SIF} Rscript ${MOTIF_R} \${vdir} ${CB} ${BIAS} ${FINEMO_WIDTH} ${ISM_FLANK} \
  || echo "WARNING: motif annotation failed for \${rsid}" >&2
for h in ${PLOT_HEADS}; do
  singularity exec ${R_SIF} Rscript ${PLOT_R} \${vdir} \${h} ${LOGO_FLANK} ${TRACK_FLANK}
done

if [[ "${MERGE_PDF}" == "1" ]]; then
  pdfs=()
  for h in ${PLOT_HEADS}; do
    f=\${vdir}/\${rsid}.contrib_\${h}.pdf; if [[ -s \${f} ]]; then pdfs+=("\${f}"); fi
  done
  for ctx in ${MERGE_ORDER}; do
    f=\${vdir}/motifs/fimo/fimo_summary_\${rsid}_\${ctx}.pdf; if [[ -s \${f} ]]; then pdfs+=("\${f}"); fi
  done
  if [[ \${#pdfs[@]} -gt 0 ]]; then
    singularity exec ${R_SIF} Rscript -e 'a <- commandArgs(TRUE); qpdf::pdf_combine(a[-1], a[1])' \
      "\${vdir}/\${rsid}.report.pdf" "\${pdfs[@]}" && echo "merged \${#pdfs[@]} PDFs"
  fi
fi
echo "done \${rsid}"
EOF

if [[ ${DRYRUN} -eq 1 ]]; then
  echo "[dry run] wrote ${splitF} and ${arrayF} (${n_found} variants)"; exit 0
fi

dep=$(printf ' && %s' "${job_deps[@]}"); dep=${dep:4}
sub=$(eval bsub ${SPLIT_BSUB} -J "snpL_${BATCH_NAME}_split" -w "'${dep}'" \
      -oo "${BATCH}/logs/split.out" -eo "${BATCH}/logs/split.err" bash "${splitF}")
echo "${sub}"
split_jid=$(sed -n 's/^Job <\([0-9]*\)>.*/\1/p' <<< "${sub}")

eval bsub ${ARRAY_BSUB} -J "\"snpL_${BATCH_NAME}[1-${n_found}]%${ARRAY_CONCURRENCY}\"" \
  -w "'done(${split_jid})'" \
  -oo "${BATCH}/logs/variant_%I.out" -eo "${BATCH}/logs/variant_%I.err" bash "${arrayF}"

echo "Submitted: 6 model jobs -> split -> ${n_found} per-variant jobs"
echo "  batch dir : ${BATCH}"
echo "  results   : ${OUTROOT}/<rsID>/"

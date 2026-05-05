#!/bin/bash
# run_snp_pipeline.sh
# Usage: bash run_snp_pipeline.sh <snp_id> [snp_bed_override]
set -euo pipefail

snp_id="${1:?Usage: run_snp_pipeline.sh <snp_id> [snp_bed_override]}"
snp_bed_override="${2:-}"

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/scripts" && pwd)"
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workingdir="${BASE_DIR}/${snp_id}"

snp_id=rs838146
SCRIPTS_DIR="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/"
BASE_DIR="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/"
workingdir="${BASE_DIR}/${snp_id}"

mkdir -p "${workingdir}"
echo "=== SNP pipeline: ${snp_id} ===" ; date
echo "Output directory: ${workingdir}"

# ── Step 1: compile BED (skip if bed already exists or override provided) ──
bash "${SCRIPTS_DIR}/01_compile_snp_bed.sh" "${workingdir}" "${snp_id}"
snp_bed=${workingdir}/${snp_id}.bed

# ── Step 2: ATAC overlap ────────────────────────────────────────────────────
echo -e "\n[Step 2] Overlapping with ATACseq peaks..." ; date
bash "${SCRIPTS_DIR}/02_overlap_atac.sh" "${workingdir}" "${snp_id}" "${snp_bed}"

# ── Step 3: TF ChIPseq overlap ──────────────────────────────────────────────
echo -e "\n[Step 3] Overlapping with TF ChIPseq peaks..." ; date
bash "${SCRIPTS_DIR}/03_overlap_tfbs.sh" "${workingdir}" "${snp_id}" "${snp_bed}"

# ── Step 4: FIMO motif analysis ─────────────────────────────────────────────
echo -e "\n[Step 4] Running FIMO..." ; date
bash "${SCRIPTS_DIR}/04_run_fimo.sh" "${workingdir}" "${snp_id}" "${snp_bed}"
ct=FRB
cond=TNF
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTS_DIR}/04_run_fimo_analyze.R "${workingdir}" "${snp_id}" "${ct}" "${cond}"

singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
        Rscript ${SCRIPTS_DIR}/04_annotate_SNP_context.R "${workingdir}" "${snp_id}"

singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript ${SCRIPTS_DIR}/04_compile_fimo_output_plots.R "${workingdir}" "${snp_id}" "${ct}" "${cond}"
echo -e "\n=== Pipeline complete for ${snp_id} ===" ; date

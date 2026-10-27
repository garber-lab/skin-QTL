#!/usr/bin/env bash
# Personal chromBPNet prediction for YOUR OWN genotyped donors: SUBMIT script.
# Same pipeline as multipeak_submit.sh (1000G), but the genotypes come from your donors'
# phased VCF (run phase_donors_beagle.sh first), and peaks can come from a list (e.g. eQTL peaks).
# Output matrices are peaks x YOUR donors, ready to compare with their measured ATAC and expression.
#
# Run on the login node (NOT with bsub):
#     bash personal_donors_submit.sh
#
# Needs multipeak_chunk.sh and multipeak_merge.sh in the same folder as this script.
#
# What it does:
# 1. One row per ATAC peak (strongest MACS2 summit), autosomes chr1-22 only
# 2. Add the summit-centered 2114 bp model window; drop windows that run off the chromosome
# 3. Pick peaks: PEAK_MODE=list (peaks overlapping PEAK_LIST), random (N_PEAKS, fixed SEED) or all
#    Also: check the donor VCF is phased, and that SAMPLES_FILE / ID_MAP match it
# 4. Split them into chunks of CHUNK_SIZE peaks
# 5. Write config.sh, copy the worker scripts into the run folder
# 6. Submit an LSF job array (one job per chunk, multipeak_chunk.sh)
#    plus a merge job (multipeak_merge.sh) that starts after all chunks have ended
#
# Each chunk job starts the container and loads each of the 10 models ONCE,
# then predicts all unique haplotypes of all its peaks.

set -eo pipefail

BASE=/pi/manuel.garber-umw/human/skin/eQTLs
CHROMBPNET=${BASE}/chromatin/chrombpnet

# ---------------- settings ----------------
RUN_NAME=KRT_donors_eQTLpeaks              # output folder name under personalGenome/multipeak/
CELLTYPE=KRT
CONDITIONS="PBS IFNG"                      # output labels
MODEL_CONDS="PBS IFN"                      # matching labels in model directory names (same order)
FOLDS="0 1 2 3 4"
PEAK_COND=IFN                              # peak file defining peaks + windows, same for every model
PEAK_MODE=list                             # list = peaks overlapping PEAK_LIST; random = N_PEAKS random; all
PEAK_LIST=/path/to/eQTL_regions.bed        # PEAK_MODE=list: BED of regions (e.g. eQTL lead SNPs + LD proxies)
N_PEAKS=100                                # PEAK_MODE=random: number of random peaks
SEED=42                                    # random seed for the peak selection
CHUNK_SIZE=25                              # peaks per LSF job
INPUT_LEN=2114
BATCH_SIZE=128                             # haplotypes per model batch (fixed size)
SAVE_FOLDS=1                               # 1 = per-fold matrices too; 0 = only fold-averaged
CLEANUP=1                                  # 1 = delete per-chunk intermediate files after success
USE_GPU=0                                  # 1 = run chunks on the gpu queue with singularity --nv

# parallelism: one LSF job = one chunk of peaks x one group of models
MODELS_PER_JOB=1                           # 1 = each model in its own job (most parallel); 10 = all models in one job
MAX_RUNNING=200                            # max jobs of each array running at once (be fair to the cluster)
MAX_ARRAY_SIZE=1000                        # your LSF MAX_JOB_ARRAY_SIZE (check: bparams -a | grep MAX_JOB_ARRAY_SIZE)

# LSF resources per job
QUEUE=short
WALL=2:00
NCPU=4
MEM=4000                                   # MB per core

PEAKS=${BASE}/ATACseq_DolphinNext/encode_atacseq_pipeline/MACS2_for_chrombpnet/ATAC_${CELLTYPE}_${PEAK_COND}_merged_alldonors_peaks.narrowPeak
FASTA=${CHROMBPNET}/ATAC_FRB_PBS_merged_alldonors/biasFactor0.5/data/downloads/hg38.fa
# your donors' PHASED genotypes (output of phase_donors_beagle.sh); @CHROM@ = chromosome name
VCF_TEMPLATE=${BASE}/chromatin/personalGenome/donor_genotypes_phased/donors_phased_@CHROM@.vcf.gz
VCF_CHROM_STYLE=chr                        # "chr" if the VCF writes chr1, "nochr" if it writes 1
SAMPLES_FILE=""                            # optional: file with the VCF sample IDs to keep (one per line); "" = all
ID_MAP=""                                  # optional: 2 columns, VCF sample ID -> donor name used in your ATAC/RNA tables
STRICT_REF=1                               # 1 = stop on REF/hg38 mismatch; 0 = skip and count such sites
SIF=${CHROMBPNET}/chrombpnet_latest.sif
RUNDIR=${BASE}/chromatin/personalGenome/multipeak/${RUN_NAME}
# -------------------------------------------

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
die() { echo "ERROR $*" >&2; exit 1; }

[[ ! -e ${RUNDIR}/config.sh ]] || die "${RUNDIR} already exists; pick a new RUN_NAME or delete it"
[[ ${PEAK_MODE} != list || -e ${PEAK_LIST} ]] || die "Missing PEAK_LIST: ${PEAK_LIST}"
for f in "${PEAKS}" "${FASTA}" "${FASTA}.fai" "${SIF}" \
         "${SCRIPT_DIR}/multipeak_chunk.sh" "${SCRIPT_DIR}/multipeak_merge.sh" "${SCRIPT_DIR}/multipeak_worker.py"; do
    [[ -e ${f} ]] || die "Missing file: ${f}"
done

# models as LABEL=path, e.g. PBS_fold0=/pi/.../ATAC_KRT_PBS_merged_alldonors/.../fold0/.../chrombpnet_nobias.h5
read -r -a CONDS <<< "${CONDITIONS}"
read -r -a MCONDS <<< "${MODEL_CONDS}"
[[ ${#CONDS[@]} -eq ${#MCONDS[@]} ]] || die "CONDITIONS and MODEL_CONDS must have the same length"
MODELS=""
for i in "${!CONDS[@]}"; do
    for F in ${FOLDS}; do
        M=${CHROMBPNET}/ATAC_${CELLTYPE}_${MCONDS[$i]}_merged_alldonors/biasFactor0.5/fold${F}/chrombpnet_model/models/chrombpnet_nobias.h5
        [[ -e ${M} ]] || die "Missing model: ${M}"
        MODELS="${MODELS} ${CONDS[$i]}_fold${F}=${M}"
    done
done

mkdir -p "${RUNDIR}/chunks" "${RUNDIR}/logs" "${RUNDIR}/merged" "${RUNDIR}/scripts"
cd "${RUNDIR}"

# ---------------- 1. one row per summit, chr1-22 ----------------
# The master-summit peak files are already de-duplicated (summits within 250 bp merged), and one wide
# peak can legitimately keep two summits > 250 bp apart (same start/end, different summit offset).
# So only rows with the SAME summit position (chrom + start + offset) are collapsed here, keeping the
# strongest: highest -log10 q (col 9), then fold enrichment (col 7).
awk -F'\t' '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' "${PEAKS}" \
    | sort -t$'\t' -k1,1 -k2,2n -k3,3n -k9,9gr -k7,7gr \
    | awk -F'\t' '!seen[$1 FS ($2 + $10)]++' > peaks_one_summit.narrowPeak
echo "Peaks (unique summits, chr1-22): $(wc -l < peaks_one_summit.narrowPeak)"

# ---------------- 2. add center + model window, drop windows off the chromosome edge ----------------
# appended columns: 11 = center (summit), 12 = window start (0-based), 13 = window end (0-based, excluded)
awk -v L="${INPUT_LEN}" 'BEGIN{FS = OFS = "\t"}
    NR == FNR { len[$1] = $2; next }
    { c = $2 + $10; ws = c - int(L / 2); we = ws + L
      # keep only the 10 narrowPeak columns (merged files may carry extra columns, e.g. both/IFN/PBS)
      if (ws >= 0 && we <= len[$1]) print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, c, ws, we }' "${FASTA}.fai" peaks_one_summit.narrowPeak > peaks_windows.tsv
echo "Peaks with a full ${INPUT_LEN} bp window: $(wc -l < peaks_windows.tsv)"

# ---------------- 3. peak selection, sorted by position ----------------
if [[ ${PEAK_MODE} == list ]]; then
    module load bedtools
    # peaks whose ATAC peak (cols 1-3) overlaps any region in PEAK_LIST
    cut -f1-3 "${PEAK_LIST}" | grep -v '^#' > peak_list_regions.bed
    bedtools intersect -u -a peaks_windows.tsv -b peak_list_regions.bed \
        | sort -t$'\t' -k1,1V -k2,2n > selected_peaks.tsv
elif [[ ${PEAK_MODE} == random ]]; then
    awk -v seed="${SEED}" 'BEGIN{srand(seed); OFS = "\t"} {print rand(), $0}' peaks_windows.tsv \
        | sort -t$'\t' -k1,1g | awk -v n="${N_PEAKS}" 'NR <= n' | cut -f2- \
        | sort -t$'\t' -k1,1V -k2,2n > selected_peaks.tsv
else
    sort -t$'\t' -k1,1V -k2,2n peaks_windows.tsv > selected_peaks.tsv
fi
N_SELECTED=$(wc -l < selected_peaks.tsv)
echo "Selected peaks (${PEAK_MODE}): ${N_SELECTED}"
[[ ${N_SELECTED} -gt 0 ]] || die "No peaks selected"

# ---------------- 3b. donor VCF checks: samples, ID map, phasing ----------------
module load bcftools
C1=$(awk 'NR == 1 {print $1}' selected_peaks.tsv)
V1=${C1}; if [[ ${VCF_CHROM_STYLE} == nochr ]]; then V1=${C1#chr}; fi
VCF1=${VCF_TEMPLATE//@CHROM@/${V1}}
[[ -e ${VCF1} ]] || die "Missing donor VCF: ${VCF1}"
bcftools query -l "${VCF1}" > vcf_samples.txt
if [[ -n ${SAMPLES_FILE} ]]; then
    MISSING=$(grep -vxF -f vcf_samples.txt "${SAMPLES_FILE}" | wc -l || true)
    [[ ${MISSING} -eq 0 ]] || die "${MISSING} SAMPLES_FILE IDs not in ${VCF1}"
fi
if [[ -n ${ID_MAP} ]]; then
    UNMAPPED=$(cut -f1 "${ID_MAP}" | grep -vxF -f - vcf_samples.txt | wc -l || true)
    echo "VCF samples without an ID_MAP entry: ${UNMAPPED} (must be 0 for the samples you keep)"
fi
echo "Donors in VCF: $(wc -l < vcf_samples.txt)"
# heterozygous genotypes on the first chromosome: phased (0|1) vs unphased (0/1)
bcftools query -f '[%GT\n]' "${VCF1}" | awk '
    { n = split($1, a, /[|\/]/); if (n == 2 && a[1] != a[2] && a[1] != "." && a[2] != ".") { het++; if ($1 ~ /\//) un++ } }
    END { printf "Heterozygous genotypes on first chromosome: %d, unphased: %d\n", het, un + 0
          if (un > 0) print "WARNING: donor VCF is not (fully) phased -> run phase_donors_beagle.sh first; results are only exact for donors with <=1 het per window (see n_donors_phase_ambiguous)" }'

# ---------------- 4. split into chunks ----------------
awk -v cs="${CHUNK_SIZE}" '{ print > ("chunks/chunk_" (int((NR - 1) / cs) + 1) ".tsv") }' selected_peaks.tsv
N_CHUNKS=$(( (N_SELECTED + CHUNK_SIZE - 1) / CHUNK_SIZE ))
echo "Chunks: ${N_CHUNKS} x up to ${CHUNK_SIZE} peaks"
N_MODELS=$(echo ${MODELS} | wc -w)
N_GROUPS=$(( (N_MODELS + MODELS_PER_JOB - 1) / MODELS_PER_JOB ))
N_JOBS=$(( N_CHUNKS * N_GROUPS ))
echo "Jobs: ${N_JOBS} (${N_CHUNKS} chunks x ${N_GROUPS} model groups of up to ${MODELS_PER_JOB} models)"

# ---------------- 5. config + copies of the worker scripts ----------------
cat > config.sh <<EOF
# written by multipeak_submit.sh on $(date '+%F %T')
RUN_NAME="${RUN_NAME}"
CELLTYPE="${CELLTYPE}"
CONDITIONS="${CONDITIONS}"
FOLDS="${FOLDS}"
PEAK_COND="${PEAK_COND}"
N_PEAKS="${N_PEAKS}"
SEED="${SEED}"
CHUNK_SIZE="${CHUNK_SIZE}"
N_CHUNKS="${N_CHUNKS}"
INPUT_LEN="${INPUT_LEN}"
VCF_CHROM_STYLE="${VCF_CHROM_STYLE}"
SAMPLES_FILE="${SAMPLES_FILE}"
ID_MAP="${ID_MAP}"
STRICT_REF="${STRICT_REF}"
PEAK_MODE="${PEAK_MODE}"
PEAK_LIST="${PEAK_LIST}"
BATCH_SIZE="${BATCH_SIZE}"
SAVE_FOLDS="${SAVE_FOLDS}"
CLEANUP="${CLEANUP}"
USE_GPU="${USE_GPU}"
PEAKS="${PEAKS}"
FASTA="${FASTA}"
VCF_TEMPLATE="${VCF_TEMPLATE}"
SIF="${SIF}"
RUNDIR="${RUNDIR}"
MODELS="${MODELS# }"
MODELS_PER_JOB="${MODELS_PER_JOB}"
N_GROUPS="${N_GROUPS}"
QUEUE="${QUEUE}"
WALL="${WALL}"
NCPU="${NCPU}"
MEM="${MEM}"
EOF
cp "${SCRIPT_DIR}/multipeak_chunk.sh" "${SCRIPT_DIR}/multipeak_merge.sh" "${SCRIPT_DIR}/multipeak_worker.py" scripts/

# ---------------- 6. submit job arrays (blocks of MAX_ARRAY_SIZE) + merge job ----------------
# job number J = OFFSET + array index; multipeak_chunk.sh turns J into (chunk, model group)
JOB=mp_${RUN_NAME}
if [[ ${USE_GPU} == 1 ]]; then QUEUE=gpu; fi
DEPS=""
B=0
for (( START = 1; START <= N_JOBS; START += MAX_ARRAY_SIZE )); do
    B=$((B + 1))
    END=$(( START + MAX_ARRAY_SIZE - 1 )); if (( END > N_JOBS )); then END=${N_JOBS}; fi
    bsub -J "${JOB}_b${B}[1-$(( END - START + 1 ))]%${MAX_RUNNING}" -q "${QUEUE}" -W "${WALL}" -n "${NCPU}" \
         -R "span[hosts=1]" -R "rusage[mem=${MEM}]" \
         -o "${RUNDIR}/logs/job_b${B}_%I.%J.log" -e "${RUNDIR}/logs/job_b${B}_%I.%J.err" \
         bash "${RUNDIR}/scripts/multipeak_chunk.sh" "${RUNDIR}" "$(( START - 1 ))"
    DEPS="${DEPS:+${DEPS} && }ended(${JOB}_b${B})"
done
bsub -J "${JOB}_merge" -w "${DEPS}" -q short -W 2:00 -n 1 -R "rusage[mem=4000]" \
     -o "${RUNDIR}/logs/merge.%J.log" -e "${RUNDIR}/logs/merge.%J.err" \
     bash "${RUNDIR}/scripts/multipeak_merge.sh" "${RUNDIR}"

echo "Submitted ${N_JOBS} jobs in ${B} array(s) (${JOB}_b1..b${B}, max ${MAX_RUNNING} running per array) and ${JOB}_merge."
echo "Outputs: ${RUNDIR}/merged/"


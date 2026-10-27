#!/usr/bin/env bash
# Multi-peak personal chromBPNet prediction: worker for ONE job = one chunk of peaks x one group of models.
# Submitted as LSF job array(s) by multipeak_submit.sh / personal_donors_submit.sh; do not run by hand
# except for debugging:
#     LSB_JOBINDEX=1 bash multipeak_chunk.sh <RUNDIR> [OFFSET]
#
# Global job number J = OFFSET + LSB_JOBINDEX  ->  chunk I = (J-1) / N_GROUPS + 1,  model group G = (J-1) % N_GROUPS + 1
#
# steps:
# 1. REF sequence of every peak's 2114 bp window in the chunk (one samtools call)
# 2. variants in all windows, all samples (one bcftools call per chromosome)
# 3. multipeak_worker.py predict: build donor haplotypes, predict the unique haplotypes with this
#    group's models -> chunks/chunk_I/pred/<model>.npy, then mark done_g<G>
# 4. the LAST group of a chunk to finish runs multipeak_worker.py finalize (fold means + tables)
#    and marks the chunk DONE. A lock directory makes sure only one job finalizes.
# Finished groups/chunks are skipped, so failed jobs can simply be resubmitted.

set -eo pipefail

RUNDIR=$1
OFFSET=${2:-0}
source "${RUNDIR}/config.sh"
# optional settings (defaults = 1000G behaviour)
VCF_CHROM_STYLE=${VCF_CHROM_STYLE:-chr}   # "chr" = VCF uses chr1; "nochr" = VCF uses 1
SAMPLES_FILE=${SAMPLES_FILE:-}            # optional: one VCF sample ID per line to keep
ID_MAP=${ID_MAP:-}                        # optional: 2 columns, VCF sample ID -> output donor name
STRICT_REF=${STRICT_REF:-1}               # 1 = stop on any REF/FASTA mismatch; 0 = skip + count them

J=$(( OFFSET + ${LSB_JOBINDEX:?run as an LSF job array element (LSB_JOBINDEX not set)} ))
I=$(( (J - 1) / N_GROUPS + 1 ))           # chunk
G=$(( (J - 1) % N_GROUPS + 1 ))           # model group within the chunk
read -r -a ALL_MODELS <<< "${MODELS}"
GROUP_MODELS=("${ALL_MODELS[@]:$(( (G - 1) * MODELS_PER_JOB )):${MODELS_PER_JOB}}")
LABELS=""
for m in "${ALL_MODELS[@]}"; do LABELS="${LABELS} ${m%%=*}"; done

BED=${RUNDIR}/chunks/chunk_${I}.tsv
OUT=${RUNDIR}/chunks/chunk_${I}
WORK=${OUT}/work_g${G}
WORKER=${RUNDIR}/scripts/multipeak_worker.py
mkdir -p "${OUT}"
LOG=${OUT}/job_g${G}.log
log() { echo "$(date '+%F %T') INFO [chunk ${I} group ${G}] $*" | tee -a "${LOG}"; }
die() { echo "$(date '+%F %T') ERROR [chunk ${I} group ${G}] $*" | tee -a "${LOG}" >&2; exit 1; }

if [[ -e ${OUT}/DONE ]]; then log "chunk already DONE, skipping"; exit 0; fi
[[ -e ${BED} ]] || die "Missing chunk file ${BED}"

module load samtools
module load bcftools
module load singularity 2>/dev/null || true
set -u
NV=""; if [[ ${USE_GPU} == 1 ]]; then NV="--nv"; fi
NTHREADS=${LSB_DJOB_NUMPROC:-4}

# ---------------- 1-3. predict this group's models (skipped if already done) ----------------
if [[ ! -e ${OUT}/done_g${G} ]]; then
    : > "${LOG}"
    T0=$(date +%s)
    rm -rf "${WORK}"; mkdir -p "${WORK}"; cd "${WORK}"
    log "$(wc -l < "${BED}") peaks on $(cut -f1 "${BED}" | uniq | tr '\n' ' '); models: ${#GROUP_MODELS[@]}"

    # 1. REF sequences of all windows
    # chunk columns: 1-10 narrowPeak, 11 center, 12 window start (0-based), 13 window end (excluded)
    awk -F'\t' '{print $1 ":" $12 + 1 "-" $13}' "${BED}" > regions.txt
    samtools faidx "${FASTA}" -r regions.txt > ref_windows.fa

    # 2. variants in all windows, per chromosome
    # one line per variant: CHROM POS ID REF ALT TYPE GT(sample1) GT(sample2) ...
    : > genotypes.tsv
    SUBSET=()
    if [[ -n ${SAMPLES_FILE} ]]; then SUBSET=(-S "${SAMPLES_FILE}"); fi
    for C in $(cut -f1 "${BED}" | uniq); do
        VC=${C}                                   # chromosome name as written in the VCF
        if [[ ${VCF_CHROM_STYLE} == nochr ]]; then VC=${C#chr}; fi
        VCF=${VCF_TEMPLATE//@CHROM@/${VC}}
        [[ -e ${VCF} ]] || die "Missing VCF: ${VCF}"
        awk -F'\t' -v vc="${VC}" -v c="${C}" 'BEGIN{OFS = "\t"} $1 == c {print vc, $12 + 1, $13}' "${BED}" > "win_${C}.txt"   # 1-based
        # sample order after the optional subset (bcftools keeps the order of the -S file)
        bcftools view -h ${SUBSET[@]+"${SUBSET[@]}"} "${VCF}" | bcftools query -l - > "samples_${C}.txt"
        if [[ -e samples.txt ]]; then
            cmp -s "samples_${C}.txt" samples.txt || die "Sample order in ${VCF} differs from previous chromosome"
        else
            cp "samples_${C}.txt" samples.txt
        fi
        # %CHROM is written back as chrN so it matches the peak file
        bcftools view ${SUBSET[@]+"${SUBSET[@]}"} -R "win_${C}.txt" "${VCF}" \
            | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE[\t%GT]\n' \
            | awk -v c="${C}" 'BEGIN{FS = OFS = "\t"} {$1 = c; print}' >> genotypes.tsv
    done
    T1=$(date +%s)
    log "bcftools/samtools: $((T1 - T0)) s; variant records: $(wc -l < genotypes.tsv); samples: $(wc -l < samples.txt)"

    # 3. haplotypes + prediction with this group's models
    if ! singularity exec ${NV} "${SIF}" python "${WORKER}" predict \
            "${BED}" ref_windows.fa genotypes.tsv samples.txt "${OUT}" \
            "${INPUT_LEN}" "${BATCH_SIZE}" "${STRICT_REF}" "${ID_MAP:-none}" "${NTHREADS}" "${G}" \
            "${GROUP_MODELS[@]}" 2>&1 | tee -a "${LOG}"; then
        die "Prediction step failed"
    fi
    T2=$(date +%s)
    printf "%s\t%s\tbcftools_samtools\t%s\t0\n%s\t%s\twall_total\t%s\t%s\n" \
        "${I}" "${G}" "$((T1 - T0))" "${I}" "${G}" "$((T2 - T0))" "$(wc -l < "${BED}")" >> "${OUT}/timing_g${G}.tsv"
    cd "${OUT}"
    if [[ ${CLEANUP} == 1 ]]; then rm -rf "${WORK}"; fi
    touch "${OUT}/done_g${G}"
    log "group done in $((T2 - T0)) s"
fi

# ---------------- 4. finalize the chunk once every model group is done ----------------
N_DONE=$(find "${OUT}" -maxdepth 1 -name 'done_g*' | wc -l)
if [[ ${N_DONE} -lt ${N_GROUPS} ]]; then
    log "${N_DONE}/${N_GROUPS} model groups done; another job will finalize the chunk"
    exit 0
fi
if ! mkdir "${OUT}/finalize.lock" 2>/dev/null; then
    log "another job is finalizing this chunk"
    exit 0
fi
log "all ${N_GROUPS} model groups done; finalizing"
if ! singularity exec "${SIF}" python "${WORKER}" finalize \
        "${OUT}" "${SAVE_FOLDS}" "${CONDITIONS}" "${LABELS# }" 2>&1 | tee -a "${LOG}"; then
    rmdir "${OUT}/finalize.lock"
    die "Finalize step failed"
fi
{ printf "chunk\tgroup\tstep\tseconds\tn\n"; cat "${OUT}"/timing_g*.tsv; } > "${OUT}/timing.tsv"
if [[ ${CLEANUP} == 1 ]]; then
    rm -rf "${OUT}/pred" "${OUT}/structure.pkl" "${OUT}"/timing_g*.tsv
fi
touch "${OUT}/DONE"
rmdir "${OUT}/finalize.lock"
log "chunk DONE"


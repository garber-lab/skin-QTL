#!/usr/bin/env bash
# Step 1: personal chromBPNet predictions for our own phased donors at tested-SNP ATAC peaks.
# Run on the login node (NOT with bsub); it submits the LSF jobs itself:
#
#     bash step1_own_donors.sh pilot            # KRT, 1000 random tested-SNP peaks (measure speed first)
#     bash step1_own_donors.sh full KRT         # all tested-SNP peaks, one cell type at a time
#                                               # (FRB, MEL after KRT's jobs finish: 5,000 pending-job cap)
#
# What it does:
# 1. Check the phased donor VCFs (chr1-22) exist and are indexed (indexes missing ones)
# 1b. Build the tested-SNP BED from the PHASED VCFs (once; reused afterwards):
#     biallelic SNVs with minor allele frequency >= MIN_MAF among the donors, saved next to the old
#     all_genotyped_SNPs.bed (which came from the old unphased VCF with mismatching REF alleles)
# 2. Per cell type: tested-SNP peaks = master-summit peaks overlapping any SNP in that BED
#    (pilot: a random 1000 of them, fixed seed)
# 3. Per cell type: write a copy of personal_donors_submit.sh with this run's settings and run it
#    -> job arrays (chunk x model) + merge job, outputs in personalGenome/multipeak/<RUN_NAME>/merged/
#
# Needs in the same folder: personal_donors_submit.sh, multipeak_chunk.sh, multipeak_worker.py, multipeak_merge.sh

set -eo pipefail

MODE=${1:?usage: bash step1_own_donors.sh pilot|full [CELLTYPE ...]}
shift
if [[ $# -gt 0 ]]; then CELLTYPES="$*"; elif [[ ${MODE} == pilot ]]; then CELLTYPES="KRT"; else CELLTYPES=""; fi

BASE=/pi/manuel.garber-umw/human/skin/eQTLs

# ---------------- settings ----------------
# phased donor genotypes (step 11 of the genotyping pipeline); @CHROM@ = chr1 ... chr22
PHASED_DIR=${BASE}/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/phased
VCF_TEMPLATE=${PHASED_DIR}/@CHROM@_imputed_combined_hg38_rsID_filtered_phased.vcf.gz
ID_MAP=""                     # optional: 2 columns, VCF sample ID -> donor name used in the RNA/ATAC tables
SAMPLES_FILE=""               # optional: VCF sample IDs to keep; "" = all donors
STRICT_REF=1                  # 1 = stop on any REF/hg38 mismatch (set 0 to skip and count them)

MACS2_DIR=${BASE}/ATACseq_DolphinNext/encode_atacseq_pipeline/MACS2_for_chrombpnet
PEAKS_TEMPLATE=${MACS2_DIR}/ATAC_@CELLTYPE@_master_summits.narrowPeak    # merged PBS+IFN summits per cell type
SNP_DIR=${BASE}/DREG/2026models/data
# OLD: ${SNP_DIR}/all_genotyped_SNPs.bed came from the old unphased VCF (REF alleles do not match) -> not used
MIN_MAF=0.05                  # "tested" SNP = minor allele frequency >= MIN_MAF among the donors; match your eQTL filter
                              # (0 = every biallelic SNV in the phased VCF)
SNP_TAG=phasedMAF${MIN_MAF}   # goes into derived file names and RUN_NAME, so old-SNP-set files are never reused
TESTED_SNPS=${SNP_DIR}/all_genotyped_SNPs.phased_hg38.maf${MIN_MAF}.bed  # built in step 1b from the phased VCFs
INPUT_DIR=${BASE}/chromatin/personalGenome/step1_inputs                  # derived peak lists go here

PILOT_N=1000                  # pilot: random tested-SNP peaks
SEED=42
PILOT_CHUNK=100               # pilot: peaks per job  (1000 peaks -> 10 chunks x 10 models = 100 jobs)
FULL_CHUNK=500                # full:  peaks per job  (~7-8 min per job on CPU)
MODELS_PER_JOB=1
# cluster limits (short queue): 512 cores per user -> 128 jobs of 4 cores; 5,000 pending jobs per user;
# MAX_JOB_ARRAY_SIZE = 10000 -> one array per cell type
MAX_RUNNING=128               # jobs running at once (= 512 cores / 4 cores per job)
MAX_ARRAY_SIZE=10000          # whole cell type in one job array
MEM=1000                      # MB per core -> 4 GB per job (benchmark peak: 1.4 GB)
MAX_PENDING=5000              # per-user pending-job cap
SAVE_FOLDS=1                  # 1 = per-fold matrices too; 0 = fold means only (~5x less disk)
USE_GPU=0
# -------------------------------------------

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
die() { echo "ERROR $*" >&2; exit 1; }
for f in personal_donors_submit.sh multipeak_chunk.sh multipeak_worker.py multipeak_merge.sh; do
    [[ -e ${SCRIPT_DIR}/${f} ]] || die "Missing ${SCRIPT_DIR}/${f}"
done
[[ ${MODE} == pilot || ${MODE} == full ]] || die "MODE must be pilot or full"
if [[ ${MODE} == full ]]; then
    N_CT=$(echo ${CELLTYPES} | wc -w)
    [[ ${N_CT} -eq 1 ]] || die "full mode: give ONE cell type (e.g. bash step1_own_donors.sh full KRT); a full cell type is ~4,600 jobs and the cluster allows ${MAX_PENDING} pending jobs per user"
fi

module load bcftools
module load bedtools

# ---------------- 1. phased donor VCFs ----------------
echo "== Checking phased VCFs in ${PHASED_DIR}"
for N in $(seq 1 22); do
    V=${VCF_TEMPLATE//@CHROM@/chr${N}}
    [[ -e ${V} ]] || die "Missing phased VCF ${V} (is phasing finished?)"
    if [[ ! -e ${V}.tbi && ! -e ${V}.csi ]]; then
        echo "   indexing $(basename "${V}")"
        bcftools index -t "${V}"
    fi
done
V1=${VCF_TEMPLATE//@CHROM@/chr1}
N_DONORS=$(bcftools query -l "${V1}" | wc -l)
for N in $(seq 2 22); do
    bcftools query -l "${VCF_TEMPLATE//@CHROM@/chr${N}}" | cmp -s - <(bcftools query -l "${V1}") \
        || die "Sample list of chr${N} differs from chr1"
done
echo "   22 chromosomes OK, ${N_DONORS} donors, same sample order"

mkdir -p "${INPUT_DIR}"

# ---------------- 1b. tested-SNP BED from the phased VCFs ----------------
if [[ ! -s ${TESTED_SNPS} ]]; then
    echo "== Building ${TESTED_SNPS}"
    echo "   biallelic SNVs, minor allele frequency >= ${MIN_MAF} among the donors (6 chromosomes at a time)"
    TMPD=$(mktemp -d "${INPUT_DIR}/snpbed.XXXXXX")
    SUBSET=(); if [[ -n ${SAMPLES_FILE} ]]; then SUBSET=(-S "${SAMPLES_FILE}"); fi
    # one BED per chromosome: chrom, POS-1, POS, rsID  (BED = 0-based start, end excluded)
    build_chr() {
        local V=${VCF_TEMPLATE//@CHROM@/chr$1}
        if ! bcftools view ${SUBSET[@]+"${SUBSET[@]}"} -v snps -m2 -M2 -q "${MIN_MAF}:minor" -Ou "${V}" \
                | bcftools query -f '%CHROM\t%POS\t%ID\n' \
                | awk 'BEGIN{OFS = "\t"} {print $1, $2 - 1, $2, $3}' > "${TMPD}/chr$1.bed"; then
            touch "${TMPD}/chr$1.failed"
        fi
    }
    for N in $(seq 1 22); do
        build_chr "${N}" &
        while (( $(jobs -rp | wc -l) >= 6 )); do sleep 2; done    # at most 6 chromosomes at a time
    done
    wait
    if ls "${TMPD}"/*.failed >/dev/null 2>&1; then
        die "building the SNP BED failed for: $(cd "${TMPD}" && ls *.failed | sed 's/.failed//' | tr '\n' ' ')"
    fi
    for N in $(seq 1 22); do cat "${TMPD}/chr${N}.bed"; done > "${TMPD}/all.bed"
    mv "${TMPD}/all.bed" "${TESTED_SNPS}"           # only appears once complete
    rm -rf "${TMPD}"
fi
echo "   tested SNPs: $(wc -l < "${TESTED_SNPS}") in $(basename "${TESTED_SNPS}")"
for CT in ${CELLTYPES}; do
    echo
    echo "================ ${CT} (${MODE}) ================"
    MASTER=${PEAKS_TEMPLATE//@CELLTYPE@/${CT}}
    [[ -e ${MASTER} ]] || die "Missing ${MASTER}"

    # ---------------- 2. tested-SNP peaks ----------------
    TESTED=${INPUT_DIR}/ATAC_${CT}_master_summits.testedSNP_${SNP_TAG}.narrowPeak
    if [[ ! -s ${TESTED} ]]; then
        bedtools intersect -u -a "${MASTER}" -b "${TESTED_SNPS}" > "${TESTED}"
    fi
    echo "   tested-SNP peaks: $(wc -l < "${TESTED}") / $(wc -l < "${MASTER}")"

    if [[ ${MODE} == pilot ]]; then
        PEAKS=${INPUT_DIR}/ATAC_${CT}_master_summits.testedSNP_${SNP_TAG}.random${PILOT_N}_seed${SEED}.narrowPeak
        awk -v seed="${SEED}" 'BEGIN{srand(seed); OFS = "\t"} {print rand(), $0}' "${TESTED}" \
            | sort -t$'\t' -k1,1g | awk -v n="${PILOT_N}" 'NR <= n' | cut -f2- > "${PEAKS}"
        RUN_NAME=${CT}_donors_${SNP_TAG}_pilot${PILOT_N}
        CHUNK=${PILOT_CHUNK}
    else
        PEAKS=${TESTED}
        RUN_NAME=${CT}_donors_${SNP_TAG}
        CHUNK=${FULL_CHUNK}
    fi
    echo "   peaks for this run: $(wc -l < "${PEAKS}")  ->  RUN_NAME=${RUN_NAME}"

    # ---------------- 3. per-run copy of personal_donors_submit.sh with these settings ----------------
    # only the FIRST assignment of each setting (the settings block) is replaced
    RUN_SCRIPT=${SCRIPT_DIR}/.step1_${RUN_NAME}.sh
    awk -v OFS="" \
        -v s_RUN_NAME="${RUN_NAME}" -v s_CELLTYPE="${CT}" -v s_PEAK_MODE="all" -v s_PEAKS="${PEAKS}" \
        -v s_CHUNK_SIZE="${CHUNK}" -v s_SEED="${SEED}" -v s_VCF_TEMPLATE="${VCF_TEMPLATE}" -v s_VCF_CHROM_STYLE="chr" \
        -v s_SAMPLES_FILE="\"${SAMPLES_FILE}\"" -v s_ID_MAP="\"${ID_MAP}\"" -v s_STRICT_REF="${STRICT_REF}" \
        -v s_MODELS_PER_JOB="${MODELS_PER_JOB}" -v s_MAX_RUNNING="${MAX_RUNNING}" -v s_SAVE_FOLDS="${SAVE_FOLDS}" \
        -v s_USE_GPU="${USE_GPU}" -v s_MAX_ARRAY_SIZE="${MAX_ARRAY_SIZE}" -v s_MEM="${MEM}" '
        BEGIN { split("RUN_NAME CELLTYPE PEAK_MODE PEAKS CHUNK_SIZE SEED VCF_TEMPLATE VCF_CHROM_STYLE SAMPLES_FILE ID_MAP STRICT_REF MODELS_PER_JOB MAX_RUNNING SAVE_FOLDS USE_GPU MAX_ARRAY_SIZE MEM", keys, " ")
                for (i in keys) want[keys[i]] = 1
                val["RUN_NAME"] = s_RUN_NAME; val["CELLTYPE"] = s_CELLTYPE; val["PEAK_MODE"] = s_PEAK_MODE
                val["PEAKS"] = s_PEAKS; val["CHUNK_SIZE"] = s_CHUNK_SIZE; val["SEED"] = s_SEED
                val["VCF_TEMPLATE"] = s_VCF_TEMPLATE; val["VCF_CHROM_STYLE"] = s_VCF_CHROM_STYLE
                val["SAMPLES_FILE"] = s_SAMPLES_FILE; val["ID_MAP"] = s_ID_MAP; val["STRICT_REF"] = s_STRICT_REF
                val["MODELS_PER_JOB"] = s_MODELS_PER_JOB; val["MAX_RUNNING"] = s_MAX_RUNNING
                val["SAVE_FOLDS"] = s_SAVE_FOLDS; val["USE_GPU"] = s_USE_GPU
                val["MAX_ARRAY_SIZE"] = s_MAX_ARRAY_SIZE; val["MEM"] = s_MEM }
        match($0, /^[A-Z_]+=/) {
            k = substr($0, 1, RLENGTH - 1)
            if ((k in want) && !(k in done)) { print k, "=", val[k], "    # set by step1_own_donors.sh"; done[k] = 1; next }
        }
        { print }
        END { for (k in want) if (!(k in done)) { print "MISSING SETTING " k > "/dev/stderr"; exit 1 } }
        ' "${SCRIPT_DIR}/personal_donors_submit.sh" > "${RUN_SCRIPT}"
    bash -n "${RUN_SCRIPT}" || die "generated ${RUN_SCRIPT} has a syntax error"

    # pending-job cap: jobs already queued + this run must stay under MAX_PENDING
    N_QUEUED=$(bjobs -u "${USER}" -noheader 2>/dev/null | wc -l || true)
    N_NEW=$(( ($(wc -l < "${PEAKS}") + CHUNK - 1) / CHUNK * 10 / MODELS_PER_JOB + 1 ))
    if (( N_QUEUED + N_NEW > MAX_PENDING )); then
        die "${N_QUEUED} jobs already queued + ~${N_NEW} new > ${MAX_PENDING} per-user cap; wait for running jobs to finish"
    fi
    echo "   jobs: ~${N_NEW} new (+${N_QUEUED} already queued), max ${MAX_RUNNING} running"

    bash "${RUN_SCRIPT}"
done

echo
echo "All submitted. Follow progress with: bjobs -A   (outputs: ${BASE}/chromatin/personalGenome/multipeak/<RUN_NAME>/merged/)"


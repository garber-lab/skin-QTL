#!/bin/bash
# Step 1b – allele-specific accessibility (ASA) validation in the 6 training donors.
# Run on the login node (NOT with bsub); it submits the LSF jobs itself:
#
#   bash asa_step1b.sh setup     # 1. find the BAMs, build the BAM table + donor ID map, copy scripts to the run folder
#   bash asa_step1b.sh submit    # 2. LSF: het sites per cell type -> allele counts per BAM -> analysis
#   bash asa_step1b.sh status    #    what is done / running
#   bash asa_step1b.sh analyze   # 3. re-run only the analysis (e.g. after changing thresholds)
#   bash asa_step1b.sh analyze --both-alleles   # sensitivity: drop donor-peaks with one-allele het sites
#                                               # (likely genotype errors) -> analysis_both_alleles/
#   bash asa_step1b.sh folds     # re-extract predictions incl. each fold (3 hets jobs; counts stay valid), then
#                                # both analyses -> per-bin concordance + held-out-chromosome check
#
# BAMs: chrM_removed_ATAC_*.nodup.bam where present, else ATAC_*.nodup.bam (KRT_IFNG); see setup step 1.
# Adding a context later (e.g. KRT_IFNG once its BAMs exist): run `setup` again, then `submit`.
# Finished parts are skipped (DONE markers), so only the new BAMs are counted, then the analysis re-runs on all.
#
# What each job does:
#   hets   (1 per cell type): asa_hets.py – each donor's het SNVs in the central 1 kb of the Step 1 peaks, taken
#          from the haplotypes the model saw, + predicted hap1/hap2 counts (PBS and IFNG fold means and each fold)
#   count  (1 per BAM):       samtools mpileup (MAPQ >= 30, BQ >= 20, no duplicates) | asa_count.py
#                             -> reads per allele per site, read pairs per haplotype per peak (each pair once)
#   analyze (1):              asa_analyze.py in the chrombpnet container -> tables + figures in analysis/
set -euo pipefail

STAGE=${1:-}

# ---------------- configuration ----------------
BASE=/pi/manuel.garber-umw/human/skin/eQTLs
BAM_DIR=${BAM_DIR:-${BASE}/ATACseq_DolphinNext/encode_atacseq_pipeline/bam_dedupped}
STEP1_TEMPLATE=${STEP1_TEMPLATE:-${BASE}/chromatin/personalGenome/multipeak/@CT@_donors_phasedMAF0.05/merged}
RUN_DIR=${RUN_DIR:-${BASE}/chromatin/personalGenome/asa_training_donors}
SIF=${SIF:-${BASE}/chromatin/chrombpnet/chrombpnet_latest.sif}
DONORS=${DONORS:-"F108 F110 F111 F25 F49 F55"}     # the 6 training donors (BAM names)
ID_MAP=${ID_MAP:-}                                 # optional: 2 columns, BAM donor -> VCF sample ID; "" = guess
PY=${PY:-"singularity exec ${SIF} python"}         # python with numpy/pandas/scipy/matplotlib
RUN_LOCAL=${RUN_LOCAL:-0}                          # 1 = run jobs directly instead of bsub (for testing)
QUEUE=short
MAPQ=30
BASEQ=20

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
S=${RUN_DIR}/scripts

if [[ ${RUN_LOCAL} != 1 ]]; then
    module load samtools
    module load singularity 2>/dev/null || true
fi

# write the job command to jobs/<name>.sh, then submit it (or run it now with RUN_LOCAL=1); prints the job id
submit_job() {   # name, dependency ("" = none), wall, mem MB, command
    local name=$1 dep=$2 wall=$3 mem=$4 cmd=$5
    mkdir -p "${RUN_DIR}/jobs"
    printf '#!/bin/bash\nset -euo pipefail\n%s\n' "${cmd}" > "${RUN_DIR}/jobs/${name}.sh"
    if [[ ${RUN_LOCAL} == 1 ]]; then
        bash "${RUN_DIR}/jobs/${name}.sh" > "${RUN_DIR}/logs/${name}.out" 2>&1 \
            || { echo "FAILED: ${name} (see logs/${name}.out)" >&2; exit 1; }
        echo local
        return
    fi
    local W=(); [[ -n ${dep} ]] && W=(-w "${dep}")
    bsub -J "${name}" ${W[@]+"${W[@]}"} -q "${QUEUE}" -W "${wall}" -n 1 -R "rusage[mem=${mem}]" \
         -o "${RUN_DIR}/logs/${name}.%J.out" -e "${RUN_DIR}/logs/${name}.%J.err" bash "${RUN_DIR}/jobs/${name}.sh" \
        | sed -n 's/.*<\([0-9]*\)>.*/\1/p'
}

case ${STAGE} in
# =====================================================================================================
setup)
    mkdir -p "${RUN_DIR}"/{scripts,logs,hets,counts,analysis}
    cp "${SCRIPT_DIR}"/asa_step1b.sh "${SCRIPT_DIR}"/asa_hets.py "${SCRIPT_DIR}"/asa_count.py "${SCRIPT_DIR}"/asa_analyze.py "${S}/"
    cd "${RUN_DIR}"

    # ---- 1. BAM table: context (= donor_CT_COND), donor, cell type, condition (IFN -> IFNG), path ----
    # file names: [chrM_removed_]ATAC_<donor>_<CT>_<PBS|IFN|IFNG>_<S1|merged>.nodup.bam
    # Per context the chrM_removed_ BAM is used if it is non-empty, else the direct .nodup.bam (KRT_IFNG, re-deduped
    # 2026-10-26). chrM reads do not matter: only autosomal het sites are counted. Empty (0-byte) files are skipped.
    printf "context\tdonor\tcell_type\tcondition\tbam\n" > bam_table.tsv
    for B in "${BAM_DIR}"/chrM_removed_ATAC_*.nodup.bam "${BAM_DIR}"/ATAC_*.nodup.bam; do
        [[ -s ${B} ]] || continue
        F=$(basename "${B}")
        if [[ ${F} =~ ^(chrM_removed_)?ATAC_([^_]+)_(KRT|FRB|MEL)_(PBS|IFNG|IFN)_.*\.nodup\.bam$ ]]; then
            D=${BASH_REMATCH[2]}; CT=${BASH_REMATCH[3]}; C=${BASH_REMATCH[4]}; [[ ${C} == IFN ]] && C=IFNG
            [[ " ${DONORS} " == *" ${D} "* ]] || continue
            cut -f1 bam_table.tsv | grep -qx "${D}_${CT}_${C}" && continue      # chrM_removed_ BAM already chosen
            printf "%s\t%s\t%s\t%s\t%s\n" "${D}_${CT}_${C}" "${D}" "${CT}" "${C}" "${B}" >> bam_table.tsv
        fi
    done
    echo "== BAMs found per context (cell type x condition), and which file type was used:"
    awk -F'\t' 'NR > 1 { k = $3 "_" $4; n[k]++; t = ($5 ~ /\/chrM_removed_/) ? "chrM_removed" : "direct nodup"; s[k] = (k in s && s[k] != t) ? "mixed" : t }
                END { for (k in n) print "   " k ": " n[k] " donors (" s[k] ")" }' bam_table.tsv | sort
    for CT in KRT FRB MEL; do for C in PBS IFNG; do
        grep -q $'\t'"${CT}"$'\t'"${C}"$'\t' bam_table.tsv || echo "   ${CT}_${C}: MISSING (re-run setup + submit once its BAMs exist)"
    done; done
    NOIDX=$(awk -F'\t' 'NR > 1 { print $5 }' bam_table.tsv | while read -r B; do
        [[ -e ${B}.bai || -e ${B%.bam}.bai ]] || echo x; done | wc -l)
    echo "   BAMs: $(($(wc -l < bam_table.tsv) - 1)); without .bai index: ${NOIDX// /} (fine: mpileup -l reads the whole BAM)"

    # ---- 2. donor ID map (BAM donor -> VCF sample ID, as in the Step 1 column headers) ----
    CT1=$(awk -F'\t' 'NR == 2 { print $3 }' bam_table.tsv)
    HDR=$( (gzip -dc "${STEP1_TEMPLATE//@CT@/${CT1}}/donor_hap_ids.tsv.gz" || true) | head -1 | tr '\t' '\n' | tail -n +2)
    if [[ -n ${ID_MAP} ]]; then
        cp "${ID_MAP}" id_map.tsv
    else
        # guess: a VCF ID equal to the donor name, or containing it as a whole word (split on _ - .)
        : > id_map.tsv
        for D in ${DONORS}; do
            M=$(echo "${HDR}" | awk -v d="${D}" '{ n = split($0, a, /[_.-]/); for (i = 1; i <= n; i++) if (a[i] == d) { print; break } }')
            if [[ $(echo "${M}" | grep -c .) == 1 ]]; then
                printf "%s\t%s\n" "${D}" "${M}" >> id_map.tsv
            else
                echo "ERROR: no unique VCF sample ID for donor ${D} (matches: ${M:-none})."
                echo "  Write a 2-column file (BAM donor, VCF sample ID) and run: ID_MAP=<file> bash $0 setup"
                echo "  VCF sample IDs (first 20): $(echo "${HDR}" | head -20 | tr '\n' ' ')"
                exit 1
            fi
        done
    fi
    while read -r D V; do
        echo "${HDR}" | grep -qx "${V}" || { echo "ERROR: ${V} (donor ${D}) is not among the Step 1 donors"; exit 1; }
    done < id_map.tsv
    echo "== Donor ID map (BAM donor -> VCF sample ID); all found among the Step 1 donors. CHECK it:"
    sed 's/^/   /' id_map.tsv
    echo "Next: bash ${S}/asa_step1b.sh submit"
    ;;
# =====================================================================================================
submit)
    cd "${RUN_DIR}"
    [[ -s bam_table.tsv && -s id_map.tsv ]] || { echo "run setup first"; exit 1; }
    for CT in $(awk -F'\t' 'NR > 1 { print $3 }' bam_table.tsv | sort -u); do
        if [[ -e hets/${CT}/DONE ]]; then echo "hets ${CT}: done"; continue; fi
        CMD="${PY} ${S}/asa_hets.py ${STEP1_TEMPLATE//@CT@/${CT}} ${CT} ${RUN_DIR}/id_map.tsv ${RUN_DIR}/hets/${CT}"
        J=$(submit_job "asa_hets_${CT}" "" 2:00 16000 "${CMD}")
        printf -v "HETS_JOB_${CT}" '%s' "${J}"      # remembered for the count jobs' dependency
        echo "hets ${CT}: job ${J}"
    done

    COUNT_JOBS=()
    while IFS=$'\t' read -r CTX D CT C B; do
        [[ ${CTX} == context ]] && continue
        if [[ -e counts/${CTX}.DONE ]]; then echo "count ${CTX}: done"; continue; fi
        CMD="samtools mpileup -l ${RUN_DIR}/hets/${CT}/${D}.positions.txt -q ${MAPQ} -Q ${BASEQ} -B -d 100000 --output-QNAME ${B} \
 | ${PY} ${S}/asa_count.py ${RUN_DIR}/hets/${CT}/${D}.sites.tsv ${RUN_DIR}/counts/${CTX} ${CTX} ${D} && touch ${RUN_DIR}/counts/${CTX}.DONE"
        HJ=HETS_JOB_${CT}; DEP=""
        [[ -n ${!HJ:-} && ${RUN_LOCAL} != 1 ]] && DEP="done(${!HJ})"
        J=$(submit_job "asa_count_${CTX}" "${DEP}" 4:00 8000 "${CMD}")
        COUNT_JOBS+=("${J}")
        echo "count ${CTX}: job ${J}"
    done < bam_table.tsv

    DEP=""
    if [[ ${RUN_LOCAL} != 1 && ${#COUNT_JOBS[@]} -gt 0 ]]; then
        DEP=$(printf "done(%s) && " "${COUNT_JOBS[@]}"); DEP=${DEP% && }
    fi
    J=$(submit_job "asa_analyze" "${DEP}" 2:00 16000 "${PY} ${S}/asa_analyze.py ${RUN_DIR}")
    echo "analyze: job ${J} (starts after the count jobs). Results: ${RUN_DIR}/analysis/summary.tsv"
    ;;
# =====================================================================================================
analyze)
    OPT=${2:-}                                     # "" or --both-alleles
    J=$(submit_job "asa_analyze${OPT:+_both_alleles}" "" 2:00 16000 "${PY} ${S}/asa_analyze.py ${RUN_DIR} ${OPT}")
    echo "analyze: job ${J}"
    ;;
# =====================================================================================================
folds)
    # asa_hets.py now also writes each fold's predictions to pred.tsv. The sites/positions files it rewrites are
    # identical, so the allele counts stay valid and only the predictions + analyses are redone.
    cd "${RUN_DIR}"
    [[ ${SCRIPT_DIR} != "${S}" ]] && cp "${SCRIPT_DIR}"/asa_step1b.sh "${SCRIPT_DIR}"/asa_hets.py "${SCRIPT_DIR}"/asa_analyze.py "${S}/"
    HJ=()
    for CT in $(awk -F'\t' 'NR > 1 { print $3 }' bam_table.tsv | sort -u); do
        CMD="${PY} ${S}/asa_hets.py ${STEP1_TEMPLATE//@CT@/${CT}} ${CT} ${RUN_DIR}/id_map.tsv ${RUN_DIR}/hets/${CT}"
        J=$(submit_job "asa_hets_${CT}" "" 2:00 16000 "${CMD}")
        HJ+=("${J}")
        echo "hets ${CT}: job ${J}"
    done
    DEP=""
    if [[ ${RUN_LOCAL} != 1 ]]; then
        DEP=$(printf "done(%s) && " "${HJ[@]}"); DEP=${DEP% && }
    fi
    J=$(submit_job "asa_analyze" "${DEP}" 2:00 16000 "${PY} ${S}/asa_analyze.py ${RUN_DIR}")
    echo "analyze: job ${J} -> analysis/"
    J=$(submit_job "asa_analyze_both_alleles" "${DEP}" 2:00 16000 "${PY} ${S}/asa_analyze.py ${RUN_DIR} --both-alleles")
    echo "analyze --both-alleles: job ${J} -> analysis_both_alleles/"
    ;;
# =====================================================================================================
status)
    cd "${RUN_DIR}"
    echo "hets done:   $(ls hets/*/DONE 2>/dev/null | wc -l) cell types"
    echo "counts done: $(ls counts/*.DONE 2>/dev/null | wc -l) / $(($(wc -l < bam_table.tsv) - 1)) BAMs"
    bjobs -w 2>/dev/null | grep asa_ || echo "no asa_ jobs in the queue"
    [[ -s analysis/summary.tsv ]] && column -t analysis/summary.tsv
    ;;
*)
    sed -n '2,20p' "$0"
    exit 1
    ;;
esac

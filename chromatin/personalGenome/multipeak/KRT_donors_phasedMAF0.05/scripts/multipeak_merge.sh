#!/usr/bin/env bash
# Multi-peak personal chromBPNet prediction: MERGE + benchmark summary.
# Submitted by the submit script to start after all job arrays have ended; can also be rerun by hand:
#     bash multipeak_merge.sh <RUNDIR>
#
# 1. List chunks without DONE and the jobs (chunk x model group) to resubmit -> merged/failed_jobs.txt
# 2. Concatenate the finished chunks' tables into merged/
# 3. Timing summary + rough extrapolation to 1k / 10k / 300k peaks -> merged/benchmark.txt

set -eo pipefail

RUNDIR=$1
source "${RUNDIR}/config.sh"
set -u
cd "${RUNDIR}"
MERGED=${RUNDIR}/merged
mkdir -p "${MERGED}"

# ---------------- 1. finished vs failed chunks / jobs ----------------
DONE_CHUNKS=()
FAILED_JOBS=()
for i in $(seq 1 "${N_CHUNKS}"); do
    if [[ -e chunks/chunk_${i}/DONE ]]; then DONE_CHUNKS+=("${i}"); continue; fi
    # rerun the groups that never finished; if all finished, rerunning group 1 retries the finalize step
    PENDING=()
    for g in $(seq 1 "${N_GROUPS}"); do
        [[ -e chunks/chunk_${i}/done_g${g} ]] || PENDING+=("$(( (i - 1) * N_GROUPS + g ))")
    done
    if [[ ${#PENDING[@]} -eq 0 ]]; then PENDING=("$(( (i - 1) * N_GROUPS + 1 ))"); fi
    FAILED_JOBS+=(${PENDING[@]+"${PENDING[@]}"})
done
printf "%s\n" ${FAILED_JOBS[@]+"${FAILED_JOBS[@]}"} | grep -v '^$' > "${MERGED}/failed_jobs.txt" || true
echo "Chunks done: ${#DONE_CHUNKS[@]} / ${N_CHUNKS}"
if [[ ${#FAILED_JOBS[@]} -gt 0 ]]; then
    LIST=$(IFS=,; echo "${FAILED_JOBS[*]}")
    echo "Jobs to rerun (global job numbers): ${LIST}  (logs: chunks/chunk_<i>/job_g<g>.log, logs/*.err)"
    echo "Resubmit only those with (if a number is above your LSF MAX_JOB_ARRAY_SIZE, use an offset instead):"
    echo "  bsub -J \"mp_${RUN_NAME}_retry[${LIST}]\" -q ${QUEUE} -W ${WALL} -n ${NCPU} -R \"span[hosts=1]\" -R \"rusage[mem=${MEM}]\" -o ${RUNDIR}/logs/retry_%I.%J.log -e ${RUNDIR}/logs/retry_%I.%J.err bash ${RUNDIR}/scripts/multipeak_chunk.sh ${RUNDIR} 0"
    echo "  then rerun: bash ${RUNDIR}/scripts/multipeak_merge.sh ${RUNDIR}"
fi
[[ ${#DONE_CHUNKS[@]} -gt 0 ]] || { echo "Nothing to merge"; exit 1; }

# ---------------- 2. concatenate tables (header from the first finished chunk) ----------------
FIRST=chunks/chunk_${DONE_CHUNKS[0]}
for f in peak_info.tsv timing.tsv; do
    { awk 'NR == 1' "${FIRST}/${f}"
      for i in "${DONE_CHUNKS[@]}"; do awk 'NR > 1' "chunks/chunk_${i}/${f}"; done
    } > "${MERGED}/${f}"
done
for path in "${FIRST}"/*.tsv.gz; do
    f=$(basename "${path}")
    { zcat "${path}" | awk 'NR == 1'
      for i in "${DONE_CHUNKS[@]}"; do zcat "chunks/chunk_${i}/${f}" | awk 'NR > 1'; done
    } | gzip > "${MERGED}/${f}"
done
echo "Merged tables: $(ls "${MERGED}" | tr '\n' ' ')"

# ---------------- 3. benchmark summary ----------------
# per job (chunk x group): fixed cost = TensorFlow import + model loading
# per chunk: compute = sum over its jobs of (wall - fixed) -> job-seconds per peak for ALL models
awk -F'\t' -v G="${N_GROUPS}" -v MPJ="${MODELS_PER_JOB}" -v NCPU="${NCPU}" '
    NR == 1 { next }
    { k = $1 SUBSEP $2 }
    $3 == "tf_import" || $3 ~ /^load_/ { fixed[k] += $4 }
    $3 ~ /^predict_/                   { pt += $4; pn += $5 }
    $3 == "parse_haplotypes"           { haps[$1] = $5 }
    $3 == "wall_total"                 { wall[k] = $4; peaks[$1] = $5; chunk[k] = $1 }
    END {
        for (k in wall) { nj++; F += fixed[k]; W += wall[k]; comp[chunk[k]] += wall[k] - fixed[k] }
        for (c in comp) { nc++; P += comp[c] / peaks[c]; NP += peaks[c]; H += haps[c] }
        F /= nj; P /= nc; speed = (pt > 0) ? pn / pt : 0
        printf "Chunks: %d   jobs: %d (%d model groups x %d models)   peaks: %d\n", nc, nj, G, MPJ, NP
        printf "Unique haplotypes: %d (%.1f per peak)\n", H, H / NP
        printf "Mean job wall time: %.0f s   fixed cost per job (TF import + model loading): %.0f s\n", W / nj, F
        printf "Job-seconds per peak (all models): %.2f   prediction speed: %.1f haplotypes/s per model (%d cores)\n\n", P, speed, NCPU
        printf "Rough extrapolation (%d models per job, %.1f unique haplotypes per peak, same resources):\n", MPJ, H / NP
        printf "%10s %12s %10s %15s %18s\n", "peaks", "chunk_size", "n_jobs", "minutes_per_job", "total_job_hours"
        split("1000 10000 300000", N, " ")
        split("25 100 500", CS, " ")
        for (i = 1; i <= 3; i++) for (j = 1; j <= 3; j++) {
            jobs = int((N[i] + CS[j] - 1) / CS[j]) * G; per = F + CS[j] * P / G
            printf "%10d %12d %10d %15.1f %18.1f\n", N[i], CS[j], jobs, per / 60, jobs * per / 3600
        }
        print "\nWall-clock time ~ total_job_hours / number of jobs running at once (set MAX_RUNNING / check your LSF limits)."
    }' "${MERGED}/timing.tsv" | tee "${MERGED}/benchmark.txt"


#!/usr/bin/env bash
# Multi-peak personal chromBPNet prediction: MERGE + benchmark summary.
# Submitted by multipeak_submit.sh to start after all chunk jobs have ended; can also be rerun by hand:
#     bash multipeak_merge.sh <RUNDIR>
#
# 1. List chunks without a DONE file -> merged/failed_chunks.txt (+ the command to resubmit them)
# 2. Concatenate the finished chunks' tables into merged/
# 3. Timing summary from all chunks + rough extrapolation to 1k / 10k / 300k peaks -> merged/benchmark.txt

set -eo pipefail

RUNDIR=$1
source "${RUNDIR}/config.sh"
set -u
cd "${RUNDIR}"
MERGED=${RUNDIR}/merged
mkdir -p "${MERGED}"

# ---------------- 1. finished vs failed chunks ----------------
DONE_CHUNKS=()
FAILED=()
for i in $(seq 1 "${N_CHUNKS}"); do
    if [[ -e chunks/chunk_${i}/DONE ]]; then DONE_CHUNKS+=("${i}"); else FAILED+=("${i}"); fi
done
printf "%s\n" "${FAILED[@]}" | grep -v '^$' > "${MERGED}/failed_chunks.txt" || true
echo "Chunks done: ${#DONE_CHUNKS[@]} / ${N_CHUNKS}"
if [[ ${#FAILED[@]} -gt 0 ]]; then
    LIST=$(IFS=,; echo "${FAILED[*]}")
    echo "Failed chunks: ${LIST} (see logs/chunk_<i>.*.err). Resubmit only those with:"
    echo "  bsub -J \"mp_${RUN_NAME}[${LIST}]\" -q short -W 2:00 -n 4 -R \"span[hosts=1]\" -R \"rusage[mem=4000]\" -o ${RUNDIR}/logs/chunk_%I.%J.log -e ${RUNDIR}/logs/chunk_%I.%J.err bash ${RUNDIR}/scripts/multipeak_chunk.sh ${RUNDIR}"
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
# per chunk: fixed cost = TensorFlow import + model loading (paid once per job);
#            per-peak cost = (wall time - fixed cost) / peaks in the chunk
awk -F'\t' '
    NR == 1 { next }
    { c = $1 }
    $2 == "tf_import" || $2 ~ /^load_/ { fixed[c] += $3 }
    $2 ~ /^predict_/                   { pred[c] += $3; npred[c] += $4; nmod[c]++ }
    $2 == "parse_haplotypes"           { haps[c] = $4 }
    $2 == "wall_total"                 { wall[c] = $3; peaks[c] = $4 }
    END {
        for (c in wall) { n++; F += fixed[c]; P += (wall[c] - fixed[c]) / peaks[c]; W += wall[c]
                          NP += peaks[c]; H += haps[c]; PT += pred[c]; PN += npred[c] }
        F /= n; P /= n
        printf "Chunks: %d   peaks: %d   unique haplotypes: %d (%.1f per peak)   models per chunk: %d\n", n, NP, H, H / NP, nmod[c]
        printf "Mean chunk wall time: %.0f s\n", W / n
        printf "Fixed cost per job (TF import + loading all models): %.0f s\n", F
        printf "Cost per peak (bcftools + haplotypes + %d model predictions + writing): %.2f s\n", nmod[c], P
        printf "Prediction speed: %.1f haplotypes/s per model\n\n", PN / PT
        printf "Rough extrapolation (same resources, %.1f unique haplotypes per peak):\n", H / NP
        printf "%10s %12s %12s %14s %18s\n", "peaks", "chunk_size", "n_jobs", "hours_per_job", "total_job_hours"
        split("1000 10000 300000", N, " ")
        for (i = 1; i <= 3; i++) for (cs = 250; cs <= 1000; cs *= 4) {
            jobs = int((N[i] + cs - 1) / cs); per = (F + cs * P) / 3600
            printf "%10d %12d %12d %14.2f %18.1f\n", N[i], cs, jobs, per, jobs * per
        }
        print "\nWall-clock time ~ total_job_hours / number of jobs running at once."
    }' "${MERGED}/timing.tsv" | tee "${MERGED}/benchmark.txt"

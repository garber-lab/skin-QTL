#!/usr/bin/env bash
# Multi-peak personal chromBPNet prediction: CHUNK worker (one LSF job-array element).
# Submitted by multipeak_submit.sh; do not run by hand except for debugging:
#     LSB_JOBINDEX=1 bash multipeak_chunk.sh <RUNDIR>
#
# steps (for the peaks in chunks/chunk_<i>.tsv):
# 1. REF sequence of every peak's 2114 bp window (one samtools call)
# 2. variants in all windows, all 1000G samples (one bcftools call per chromosome)
# 3. Python in the chrombpnet container: per peak, build donor haplotypes (biallelic SNVs);
#    collect the unique haplotypes of ALL peaks, then load each model once and predict them
# 4. mark DONE, remove intermediate files
# A chunk with a DONE file is skipped, so failed chunks can simply be resubmitted.

set -eo pipefail

RUNDIR=$1
source "${RUNDIR}/config.sh"
I=${LSB_JOBINDEX:?run as an LSF job array element (LSB_JOBINDEX not set)}
BED=${RUNDIR}/chunks/chunk_${I}.tsv
OUT=${RUNDIR}/chunks/chunk_${I}

mkdir -p "${OUT}"
LOG=${OUT}/chunk.log
log() { echo "$(date '+%F %T') INFO [chunk ${I}] $*" | tee -a "${LOG}"; }
die() { echo "$(date '+%F %T') ERROR [chunk ${I}] $*" | tee -a "${LOG}" >&2; exit 1; }

if [[ -e ${OUT}/DONE ]]; then log "already DONE, skipping"; exit 0; fi
: > "${LOG}"
[[ -e ${BED} ]] || die "Missing chunk file ${BED}"

module load samtools
module load bcftools
module load singularity 2>/dev/null || true
set -u
cd "${OUT}"
T0=$(date +%s)
log "$(wc -l < "${BED}") peaks on $(cut -f1 "${BED}" | uniq | tr '\n' ' ')"

# ---------------- 1. REF sequences of all windows ----------------
# chunk columns: 1-10 narrowPeak, 11 center, 12 window start (0-based), 13 window end (excluded)
awk -F'\t' '{print $1 ":" $12 + 1 "-" $13}' "${BED}" > regions.txt
samtools faidx "${FASTA}" -r regions.txt > ref_windows.fa

# ---------------- 2. variants in all windows, per chromosome ----------------
# one line per variant: CHROM POS ID REF ALT TYPE GT(sample1) GT(sample2) ...
: > genotypes.tsv
rm -f samples.txt
for C in $(cut -f1 "${BED}" | uniq); do
    VCF=${VCF_TEMPLATE//@CHROM@/${C}}
    [[ -e ${VCF} ]] || die "Missing VCF: ${VCF}"
    awk -F'\t' -v c="${C}" 'BEGIN{OFS = "\t"} $1 == c {print $1, $12 + 1, $13}' "${BED}" > "win_${C}.txt"   # 1-based
    if [[ -e samples.txt ]]; then
        bcftools query -l "${VCF}" | cmp -s - samples.txt || die "Sample order in ${VCF} differs from previous chromosome"
    else
        bcftools query -l "${VCF}" > samples.txt
    fi
    bcftools view -R "win_${C}.txt" "${VCF}" \
        | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE[\t%GT]\n' >> genotypes.tsv
done
T1=$(date +%s)
log "bcftools/samtools: $((T1 - T0)) s; variant records: $(wc -l < genotypes.tsv); samples: $(wc -l < samples.txt)"

# ---------------- 3. haplotypes + prediction with every model ----------------
NV=""; if [[ ${USE_GPU} == 1 ]]; then NV="--nv"; fi
read -r -a MODEL_ARR <<< "${MODELS}"
if ! singularity exec ${NV} "${SIF}" python - \
        "${BED}" ref_windows.fa genotypes.tsv samples.txt "${OUT}" \
        "${INPUT_LEN}" "${BATCH_SIZE}" "${SAVE_FOLDS}" "${I}" "${MODEL_ARR[@]}" \
        <<'PY' 2>&1 | tee -a "${LOG}"; then
import sys, time, gzip, bisect
t0 = time.time()
import numpy as np
import tensorflow as tf
timing = [("tf_import", time.time() - t0, 0)]

# --- 3a. arguments from bash ---
bed, ref_fa, geno_tsv, samples_txt, outdir, input_len, batch_size, save_folds, chunk_id = sys.argv[1:10]
input_len, batch_size = int(input_len), int(batch_size)
models = [m.split("=", 1) for m in sys.argv[10:]]        # [(label, path)], label = CONDITION_foldN
conditions = list(dict.fromkeys(lab.rsplit("_fold", 1)[0] for lab, _ in models))
print(f"TensorFlow {tf.__version__}, numpy {np.__version__}, GPUs: {len(tf.config.list_physical_devices('GPU'))}")
t1 = time.time()

# --- 3b. peaks, REF sequences, samples ---
peaks = []
for line in open(bed):
    f = line.rstrip("\n").split("\t")
    peaks.append(dict(chrom=f[0], start=int(f[1]), end=int(f[2]), name=f[3], summit=int(f[9]),
                      center=int(f[10]), ws=int(f[11])))
seqs, cur = [], None
for line in open(ref_fa):
    if line.startswith(">"):
        if cur is not None:
            seqs.append("".join(cur).upper())
        cur = []
    else:
        cur.append(line.strip())
seqs.append("".join(cur).upper())
if len(seqs) != len(peaks):
    sys.exit(f"ERROR {len(seqs)} REF sequences for {len(peaks)} peaks")
for p, s in zip(peaks, seqs):
    if len(s) != input_len:
        sys.exit(f"ERROR {p['name']}: REF length {len(s)} != {input_len}")
    p["seq"] = s
samples = [l.strip() for l in open(samples_txt) if l.strip()]
n_samples = len(samples)

# --- 3c. variant lines grouped by chromosome, sorted by position, de-duplicated ---
by_chrom = {}
for line in open(geno_tsv):
    chrom, pos, vid, ref, alt, rest = line.split("\t", 5)
    by_chrom.setdefault(chrom, {})[(int(pos), ref, alt)] = line     # overlapping windows -> keep one copy
for c in by_chrom:
    items = sorted(by_chrom[c].items())
    by_chrom[c] = ([k[0] for k, _ in items], [l for _, l in items])

# --- 3d. per peak: donor haplotypes from biallelic SNVs; collect unique haplotypes ---
lut = {b: i for i, b in enumerate("ACGT")}
entries = []                                    # global list of (peak index, tuple of ALT variants)
for pi, p in enumerate(peaks):
    pos_list, lines = by_chrom.get(p["chrom"], ([], []))
    lo = bisect.bisect_left(pos_list, p["ws"] + 1)
    hi = bisect.bisect_left(pos_list, p["ws"] + input_len + 1)
    variants, hap_vars = [], [[[], []] for _ in range(n_samples)]
    n_indel = n_multi = 0
    unphased = np.zeros(n_samples, bool)
    missing = np.zeros(n_samples, bool)
    for line in lines[lo:hi]:
        f = line.rstrip("\n").split("\t")
        pos, ref, alt, vtype, gts = int(f[1]), f[3].upper(), f[4].upper(), f[5], f[6:]
        if "," in alt:
            n_multi += 1
            continue
        if vtype != "SNP" or len(ref) != 1 or alt not in lut:
            n_indel += 1
            continue
        if len(gts) != n_samples:
            sys.exit(f"ERROR {len(gts)} genotypes at {p['chrom']}:{pos}, expected {n_samples}")
        off = pos - 1 - p["ws"]
        if p["seq"][off] != ref:
            sys.exit(f"ERROR REF mismatch at {p['chrom']}:{pos}: VCF {ref} vs FASTA {p['seq'][off]}")
        v = len(variants)
        variants.append((pos, ref, alt, off))
        for s, gt in enumerate(gts):
            a = gt.replace("/", "|").split("|")
            if len(a) == 1:
                a = a * 2
            if "." in a:
                missing[s] = True                # missing allele -> left as REF
            al = [x == "1" for x in a]
            if al[0] != al[1] and "|" not in gt:
                unphased[s] = True               # het without phase: order taken as written
            for h in (0, 1):
                if al[h]:
                    hap_vars[s][h].append(v)
    hap_index = {(): 0}                          # local hap_id 0 = hg38 REF
    hap_ids = np.zeros((n_samples, 2), dtype=np.int32)
    for s in range(n_samples):
        for h in (0, 1):
            key = tuple(hap_vars[s][h])
            if key not in hap_index:
                hap_index[key] = len(hap_index)
            hap_ids[s, h] = hap_index[key]
    keys = sorted(hap_index, key=hap_index.get)
    p.update(variants=variants, hap_ids=hap_ids, keys=keys, offset=len(entries),
             n_indel=n_indel, n_multi=n_multi, n_unphased=int(unphased.sum()), n_missing=int(missing.sum()),
             n_carriers=np.bincount(hap_ids.ravel(), minlength=len(keys)))
    oh = np.zeros((input_len, 4), dtype=np.float32)
    for i, b in enumerate(p["seq"]):
        if b in lut:
            oh[i, lut[b]] = 1.0
    p["ref_oh"] = oh
    entries += [(pi, k) for k in keys]
n_entries = len(entries)
timing.append(("parse_haplotypes", time.time() - t1, n_entries))
print(f"Peaks: {len(peaks)}; samples: {n_samples}; unique haplotypes (incl. REF per peak): {n_entries}")

def encode(pi, key):
    p = peaks[pi]
    x = p["ref_oh"].copy()
    for v in key:
        off, alt = p["variants"][v][3], p["variants"][v][2]
        x[off] = 0.0
        x[off, lut[alt]] = 1.0
    return x

def load(path):
    try:
        return tf.keras.models.load_model(path, compile=False)
    except Exception as e:
        print(f"Plain load failed ({e}); retrying with custom_objects")
        return tf.keras.models.load_model(path, compile=False,
                                          custom_objects={"tf": tf, "multinomial_nll": lambda *a: 0.0})

# --- 3e. predict all unique haplotypes with each model, fixed batch size (last batch padded) ---
counts = {}
for lab, path in models:
    t = time.time()
    model = load(path)
    if model.input_shape[1] != input_len:
        sys.exit(f"ERROR {lab}: model input length {model.input_shape[1]} != {input_len}")
    timing.append((f"load_{lab}", time.time() - t, 0))
    t = time.time()
    lc_all = []
    for b0 in range(0, n_entries, batch_size):
        batch = entries[b0:b0 + batch_size]
        x = np.zeros((batch_size, input_len, 4), dtype=np.float32)
        for j, (pi, key) in enumerate(batch):
            x[j] = encode(pi, key)
        out = model.predict_on_batch(x)
        lc_all.append(np.asarray(out[1]).reshape(-1)[:len(batch)])
    # ChromBPNet's count head is trained on log(1 + total reads), so expm1 converts it back to reads
    counts[lab] = np.expm1(np.concatenate(lc_all))
    dt = time.time() - t
    timing.append((f"predict_{lab}", dt, n_entries))
    print(f"  {lab}: load {timing[-2][1]:.1f} s, predict {dt:.1f} s ({n_entries / max(dt, 1e-9):.1f} haplotypes/s)")
    tf.keras.backend.clear_session()
    del model

for c in conditions:
    counts[f"{c}_mean"] = np.mean([counts[lab] for lab, _ in models if lab.rsplit("_fold", 1)[0] == c], axis=0)
cols = []
for c in conditions:
    cols += [lab for lab, _ in models if lab.rsplit("_fold", 1)[0] == c] + [f"{c}_mean"]
matrix_cols = cols if save_folds == "1" else [f"{c}_mean" for c in conditions]

# --- 3f. write outputs ---
t = time.time()
with open(f"{outdir}/peak_info.tsv", "w") as fh:
    fh.write("peak_id\tchrom\tstart\tend\tsummit_offset\tcenter\twindow\tn_snvs\tn_skipped_indel\t"
             "n_skipped_multiallelic\tn_unique_haps\tn_donors_unphased\tn_donors_missing\n")
    for p in peaks:
        fh.write(f"{p['name']}\t{p['chrom']}\t{p['start']}\t{p['end']}\t{p['summit']}\t{p['center']}\t"
                 f"{p['chrom']}:{p['ws'] + 1}-{p['ws'] + input_len}\t{len(p['variants'])}\t{p['n_indel']}\t"
                 f"{p['n_multi']}\t{len(p['keys'])}\t{p['n_unphased']}\t{p['n_missing']}\n")

# one row per unique haplotype of each peak (hap_id 0 = hg38 REF)
with gzip.open(f"{outdir}/haplotypes.tsv.gz", "wt") as fh:
    fh.write("\t".join(["peak_id", "hap_id", "n_donor_haplotypes", "n_alt", "alt_alleles"] + cols) + "\n")
    for p in peaks:
        for k, key in enumerate(p["keys"]):
            g = p["offset"] + k
            alts = ",".join(f"{p['variants'][v][0]}{p['variants'][v][1]}>{p['variants'][v][2]}" for v in key) or "REF"
            fh.write("\t".join([p["name"], str(k), str(p["n_carriers"][k]), str(len(key)), alts]
                               + [f"{counts[c][g]:.4f}" for c in cols]) + "\n")

# each donor's two hap_ids per peak ("hap1|hap2"), to link donors to haplotypes.tsv.gz
with gzip.open(f"{outdir}/donor_hap_ids.tsv.gz", "wt") as fh:
    fh.write("peak_id\t" + "\t".join(samples) + "\n")
    for p in peaks:
        fh.write(p["name"] + "\t" + "\t".join(f"{a}|{b}" for a, b in p["hap_ids"]) + "\n")

# peaks x donors matrices: diploid predicted counts = mean of the two haplotypes
for c in matrix_cols:
    with gzip.open(f"{outdir}/diploid_{c}.tsv.gz", "wt") as fh:
        fh.write("peak_id\t" + "\t".join(samples) + "\n")
        for p in peaks:
            g = p["offset"] + p["hap_ids"]
            dip = (counts[c][g[:, 0]] + counts[c][g[:, 1]]) / 2
            fh.write(p["name"] + "\t" + "\t".join(f"{x:.2f}" for x in dip) + "\n")
timing.append(("write_outputs", time.time() - t, 0))

with open(f"{outdir}/timing.tsv", "w") as fh:
    fh.write("chunk\tstep\tseconds\tn_haplotypes\n")
    for step, sec, n in timing:
        fh.write(f"{chunk_id}\t{step}\t{sec:.2f}\t{n}\n")
print(f"Wrote outputs; python total {time.time() - t0:.1f} s")
PY
    die "Prediction step failed"
fi

# ---------------- 4. finish ----------------
T2=$(date +%s)
printf "%s\tbcftools_samtools\t%s\t0\n%s\twall_total\t%s\t%s\n" \
    "${I}" "$((T1 - T0))" "${I}" "$((T2 - T0))" "$(wc -l < "${BED}")" >> timing.tsv
if [[ ${CLEANUP} == 1 ]]; then
    rm -f regions.txt ref_windows.fa genotypes.tsv win_*.txt
fi
touch DONE
log "DONE in $((T2 - T0)) s"

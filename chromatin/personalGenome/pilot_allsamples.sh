#!/usr/bin/env bash
#BSUB -J pilot_allsamples_KRT_IFNG
#BSUB -q short
#BSUB -W 2:00
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4000]"
#BSUB -o pilot_allsamples_KRT_IFNG.%J.log
#BSUB -e pilot_allsamples_KRT_IFNG.%J.err

# pilot: haplotype-based personal chromBPNet prediction for ALL 1000G donors and 1 ATAC-seq peak
#        using chromBPNet model KRT IFNG fold0
# steps:
# 1. SNP -> overlapping peak (strongest summit) -> 2114 bp model input window centered on the summit
# 2. Check the chromosome exists in the VCF, list all samples
# 3. Subset the phased 1000G VCF to the window, ALL samples (biallelic SNVs) -> genotype table
# 4. Check VCF REF alleles against the FASTA
# 5. Extract the 2114 bp REF sequence
# 6. In Python (inside the chrombpnet container):
#    build hap1/hap2 for every donor from the genotype table, predict each UNIQUE haplotype once,
#    write per-haplotype and per-donor (diploid) predicted counts
#
# Difference from the 1-sample script: no bcftools consensus loop per donor.
# Haplotypes are built in Python by swapping ALT bases into the REF sequence (SNVs only),
# and identical haplotypes shared by many donors are predicted only once.
#
# Usage:
# bsub < pilot_allsamples.sh

set -eo pipefail

BASE=/pi/manuel.garber-umw/human/skin/eQTLs
CHROMBPNET=${BASE}/chromatin/chrombpnet

# ---------------- settings ----------------
CELLTYPE=KRT
CONDITION=IFNG                        # filename-safe label written into output tables
MODEL_COND=IFN                            # label used in model/peak directory names
FOLD=0
SNP=rs6844475         # strong eQTL in a KRT IFNG peak with strong chromBPNet scores
SNP_REF=C             # dbSNP REF allele, checked against hg38 below
INPUT_LEN=2114
BATCH_SIZE=64         # haplotypes per model.predict batch
RUN_MODEL=1           # 0 = stop after building the genotype table and REF sequence

MODEL=${CHROMBPNET}/ATAC_${CELLTYPE}_${MODEL_COND}_merged_alldonors/biasFactor0.5/fold${FOLD}/chrombpnet_model/models/chrombpnet_nobias.h5
PEAKS=${BASE}/ATACseq_DolphinNext/encode_atacseq_pipeline/MACS2_for_chrombpnet/ATAC_${CELLTYPE}_${MODEL_COND}_merged_alldonors_peaks.narrowPeak
FASTA=${CHROMBPNET}/ATAC_FRB_PBS_merged_alldonors/biasFactor0.5/data/downloads/hg38.fa
METADATA=${BASE}/literature/1000G/GRCh38_30x_phased/metadata/20130606_g1k_3202_samples_ped_population.txt
ALLSNPBED=${BASE}/DREG/2026models/data/all_genotyped_SNPs.bed
SIF=${CHROMBPNET}/chrombpnet_latest.sif
OUTDIR=${BASE}/chromatin/personalGenome/pilot/${SNP}_allsamples_${CELLTYPE}_${CONDITION}_fold${FOLD}
# -------------------------------------------

mkdir -p "${OUTDIR}"
LOG=${OUTDIR}/pilot.log
: > "${LOG}"
log() { echo "$(date '+%F %T') INFO $*" | tee -a "${LOG}"; }
die() { echo "$(date '+%F %T') ERROR $*" | tee -a "${LOG}" >&2; exit 1; }

module load samtools
module load bcftools
module load bedtools
module load singularity 2>/dev/null || true      # skip if singularity is already on PATH
set -u                                    # after module load: unset variables are now errors
cd "${OUTDIR}"

for f in "${PEAKS}" "${FASTA}" "${FASTA}.fai" "${ALLSNPBED}" "${METADATA}"; do
    [[ -e ${f} ]] || die "Missing file: ${f}"
done
if [[ ${RUN_MODEL} == 1 ]]; then
    [[ -e ${MODEL} ]] || die "Missing model: ${MODEL}"
    [[ -e ${SIF} ]] || die "Missing container: ${SIF}"
    command -v singularity >/dev/null || die "singularity not found"
fi

# ---------------- 1. SNP -> peak -> summit -> model input window ----------------
# (same as the 1-sample script)

# 1a. Get the SNP's BED line (0-based start, end not included)
SNPBED=${OUTDIR}/snp.bed
grep -w "${SNP}" "${ALLSNPBED}" > "${SNPBED}" || true      # no match is handled just below
N_SNP=$(wc -l < "${SNPBED}")
[[ ${N_SNP} -eq 1 ]] || die "Expected 1 line for ${SNP} in ${ALLSNPBED}, found ${N_SNP}"

# 1b. Pull out the SNP's chr and pos
CHROM=$(cut -f1 "${SNPBED}")              # chr4
SNP_START=$(cut -f2 "${SNPBED}")          # 0-based: 25020561
SNP_POS=$((SNP_START + 1))                # 1-based (VCF/dbSNP): 25020562
log "SNP ${SNP} at ${CHROM}:${SNP_POS} (1-based)"

# check hg38 base at the SNP matches the expected REF allele
FASTA_BASE=$(samtools faidx "${FASTA}" "${CHROM}:${SNP_POS}-${SNP_POS}" | awk 'NR==2' | tr a-z A-Z)
[[ ${FASTA_BASE} == "${SNP_REF}" ]] || die "hg38 base at ${CHROM}:${SNP_POS} is ${FASTA_BASE}, expected ${SNP_REF}"

# 1c. Now that the chr is known, pick the 1000G VCF for it
VCF=${BASE}/literature/1000G/GRCh38_30x_phased/vcf/CCDG_14151_B01_GRM_WGS_2020-08-05_${CHROM}.filtered.shapeit2-duohmm-phased.vcf.gz
[[ -e ${VCF} ]] || die "Missing VCF: ${VCF}"

# 1d. Find the ATAC peak that contains the SNP (-wa = print the full original peak row)
PEAKBED_ALL=${OUTDIR}/peak_all_summits.bed
PEAKBED=${OUTDIR}/peak.bed
bedtools intersect -wa -a "${PEAKS}" -b "${SNPBED}" > "${PEAKBED_ALL}"
N_PEAK=$(wc -l < "${PEAKBED_ALL}")
[[ ${N_PEAK} -ge 1 ]] || die "${SNP} does not overlap any peak in ${PEAKS}"

# MACS2 --call-summits writes one row per summit for the same peak.
# Keep the strongest summit: highest -log10 q (col 9), then fold enrichment (col 7).
sort -t$'\t' -k9,9gr -k7,7gr "${PEAKBED_ALL}" | awk 'NR==1' > "${PEAKBED}"
log "${SNP} overlaps ${N_PEAK} peak row(s); using strongest summit: $(cut -f4,10 "${PEAKBED}" | tr '\t' ' ')"

# 1e. Read the peak columns into variables
P_CHROM=$(cut -f1 "${PEAKBED}")           # chr4
P_START=$(cut -f2 "${PEAKBED}")           # 25020188 (0-based)
P_END=$(cut -f3 "${PEAKBED}")             # 25020755
P_NAME=$(cut -f4 "${PEAKBED}")            # ATAC_KRT_IFN_merged_alldonors_peak_221823
P_SUMMIT=$(cut -f10 "${PEAKBED}")         # 227 = offset from P_START to the summit

# 1f. Summit position = center of the model window
[[ ${P_SUMMIT} -ge 0 ]] || die "Peak ${P_NAME} has no summit offset (col10=${P_SUMMIT})"
CENTER=$((P_START + P_SUMMIT))            # 25020188 + 227 = 25020415

# Save the summit as a 1 bp BED (for IGV / bedtools closest)
PEAKSUMMIT=${OUTDIR}/summit.bed
printf "%s\t%s\t%s\t%s\n" "${P_CHROM}" "${CENTER}" "$((CENTER + 1))" "${P_NAME}" > "${PEAKSUMMIT}"

# 1g. Model input window: INPUT_LEN bp centered on the summit
WIN_START=$((CENTER - INPUT_LEN / 2))     # 0-based start          -> 25019358
WIN_END=$((WIN_START + INPUT_LEN))        # 0-based end (excluded) -> 25021472
REGION="${P_CHROM}:$((WIN_START + 1))-${WIN_END}"   # 1-based, for samtools/bcftools -> chr4:25019359-25021472

# 1h. Sanity checks: window inside the chromosome, SNP inside the window
CHROM_LEN=$(awk -v c="${P_CHROM}" '$1 == c {print $2}' "${FASTA}.fai")
[[ -n ${CHROM_LEN} ]] || die "${P_CHROM} not in FASTA index"
(( WIN_START >= 0 && WIN_END <= CHROM_LEN )) || die "Input window runs off chromosome edge"
(( SNP_START >= WIN_START && SNP_START < WIN_END )) || die "SNP is outside the model input window"

SNP_DIST=$((SNP_START - CENTER))          # +146 = SNP is 146 bp downstream of summit
SNP_SEQ_POS=$((SNP_START - WIN_START + 1))   # position of SNP within the 2114 bp sequence (1-based) -> 1204

log "Peak ${P_NAME} ${P_CHROM}:${P_START}-${P_END} summit_offset=${P_SUMMIT} center=${CENTER}"
log "Input window ${REGION} (1-based), length ${INPUT_LEN}"
log "SNP is ${SNP_DIST} bp from summit; base ${SNP_SEQ_POS} of the input sequence"


# ---------------- 2. chromosome check, list all samples ----------------
bcftools view -h "${VCF}" | grep "^##contig=<ID=${P_CHROM}[,>]" >/dev/null || die "${P_CHROM} not in VCF contigs (chr1 vs 1?)"
N_VCF_SAMPLES=$(bcftools query -l "${VCF}" | wc -l)
log "Samples in VCF: ${N_VCF_SAMPLES}"


# ---------------- 3. variants in 2114bp window, ALL samples ----------------
# 3a. whole chr VCF -> only the 2114bp window, all samples
ALL_VCF=${OUTDIR}/window_all.vcf.gz
bcftools view -r "${REGION}" -Oz -o "${ALL_VCF}" "${VCF}"
bcftools index -f -t "${ALL_VCF}"

# 3b. biallelic SNVs only (site-level filter).
#     No -p / GT filters here: those act per sample. Unphased or missing genotypes
#     are counted per donor in step 6 instead.
SNV_VCF=${OUTDIR}/window_snvs.vcf.gz
bcftools view -v snps -m2 -M2 -Oz -o "${SNV_VCF}" "${ALL_VCF}"
bcftools index -f -t "${SNV_VCF}"

N_ALL=$(bcftools view -H "${ALL_VCF}" | wc -l)
N_SNV=$(bcftools view -H "${SNV_VCF}" | wc -l)
log "Variant records in window: ${N_ALL}; biallelic SNVs used: ${N_SNV}"
log "Skipped: indel=$(bcftools view -H -V snps "${ALL_VCF}" | wc -l) \
multiallelic_snv=$(bcftools view -H -v snps -m3 "${ALL_VCF}" | wc -l)"

# 3c. genotype table: one row per SNV, one GT column per sample (sample order in samples.txt)
SAMPLES_TXT=${OUTDIR}/samples.txt
GENO_TSV=${OUTDIR}/genotypes.tsv
bcftools query -l "${SNV_VCF}" > "${SAMPLES_TXT}"
bcftools query -f '%POS\t%ID\t%REF\t%ALT[\t%GT]\n' "${SNV_VCF}" > "${GENO_TSV}"
log "Genotype table: $(wc -l < "${GENO_TSV}") SNVs x $(wc -l < "${SAMPLES_TXT}") samples"

# 3d. check SNP is among the SNVs, and log its genotype counts across donors
SNP_ALT=$(bcftools query -r "${CHROM}:${SNP_POS}" -f '%ALT\n' "${SNV_VCF}" | awk 'NR==1')
[[ -n ${SNP_ALT} ]] || die "${SNP} (${CHROM}:${SNP_POS}) is not a biallelic SNV in ${VCF}"
log "${SNP} ${SNP_REF}>${SNP_ALT}; genotype counts across donors: \
$(bcftools query -r "${CHROM}:${SNP_POS}" -f '[%GT\n]' "${SNV_VCF}" | sort | uniq -c | awk '{printf "%s=%s ", $2, $1}')"


# ---------------- 4. REF check ----------------
# bcftools norm --check-ref e exits with an error on any REF/FASTA mismatch
bcftools norm -f "${FASTA}" --check-ref e -Ou "${ALL_VCF}" > /dev/null 2>"${OUTDIR}/refcheck.err" \
    || die "REF allele mismatch vs FASTA (see ${OUTDIR}/refcheck.err)"
log "REF alleles match FASTA"


# ---------------- 5. REF sequence of the window ----------------
# one plain uppercase line; donors' haplotypes are built from it in step 6
REF_TXT=${OUTDIR}/ref_seq.txt
samtools faidx "${FASTA}" "${REGION}" | grep -v '^>' | tr -d '\n' | tr 'acgtn' 'ACGTN' > "${REF_TXT}"
REF_LEN=$(wc -c < "${REF_TXT}")
[[ ${REF_LEN} -eq ${INPUT_LEN} ]] || die "REF sequence length ${REF_LEN} != ${INPUT_LEN}"
log "REF sequence: ${REF_LEN} bp"

if [[ ${RUN_MODEL} != 1 ]]; then
    log "RUN_MODEL=0; stopping before prediction. Outputs in ${OUTDIR}/"
    exit 0
fi


# ---------------- 6. build haplotypes + ChromBPNet prediction (inside the chrombpnet container) ----------------
# The Python code between <<'PY' and PY is fed straight into the container's python and run here.
# Arguments after "python -" arrive in Python as sys.argv[1:].

log "Predicting with ${MODEL}"
if ! singularity exec "${SIF}" python - \
        "${MODEL}" "${REF_TXT}" "${GENO_TSV}" "${SAMPLES_TXT}" "${METADATA}" "${OUTDIR}" \
        "${INPUT_LEN}" "${WIN_START}" "${CENTER}" "${SNP_POS}" "${BATCH_SIZE}" \
        "${P_NAME}" "${CELLTYPE}" "${CONDITION}" "${FOLD}" \
        <<'PY' 2>&1 | tee -a "${LOG}"; then
import sys
import numpy as np
import tensorflow as tf

# --- 6a. arguments from bash ---
(model_path, ref_txt, geno_tsv, samples_txt, metadata, outdir,
 input_len, win_start, center, snp_pos, batch_size,
 peak, celltype, condition, fold) = sys.argv[1:]
input_len, win_start, center, snp_pos, batch_size = map(int, (input_len, win_start, center, snp_pos, batch_size))
print(f"TensorFlow {tf.__version__}, numpy {np.__version__}")

# --- 6b. load model and check it expects INPUT_LEN bp ---
try:
    model = tf.keras.models.load_model(model_path, compile=False)
except Exception as e:
    print(f"Plain load failed ({e}); retrying with custom_objects")
    model = tf.keras.models.load_model(model_path, compile=False,
                                       custom_objects={"tf": tf, "multinomial_nll": lambda *a: 0.0})
print(f"Model input shape {model.input_shape}; outputs {[tuple(o.shape) for o in model.outputs]}")
if model.input_shape[1] != input_len:
    sys.exit(f"ERROR model input length {model.input_shape[1]} != INPUT_LEN {input_len}")

# --- 6c. read REF sequence, samples, population metadata ---
ref_seq = open(ref_txt).read().strip()
if len(ref_seq) != input_len:
    sys.exit(f"ERROR REF sequence length {len(ref_seq)} != {input_len}")
samples = [l.strip() for l in open(samples_txt) if l.strip()]
n_samples = len(samples)

pop = {}
for line in open(metadata):
    f = line.split()
    if len(f) >= 7 and f[1] != "SampleID":
        pop[f[1]] = (f[5], f[6])

# --- 6d. read genotype table -> which variants each haplotype carries ---
# a haplotype is described by the list of variant indices where it carries ALT
variants = []                                   # (pos, id, ref, alt, offset in window)
hap_vars = [[[], []] for _ in range(n_samples)] # hap_vars[sample][0 or 1] = [variant indices]
n_unphased = np.zeros(n_samples, dtype=int)     # heterozygous genotypes written with "/"
n_missing = np.zeros(n_samples, dtype=int)      # genotypes containing "."
snp_dosage = np.full(n_samples, -1)             # 0/1/2 ALT copies of the focal SNP (-1 = not found)
snp_gt = ["NA"] * n_samples

for line in open(geno_tsv):
    f = line.rstrip("\n").split("\t")
    pos, vid, ref, alt, gts = int(f[0]), f[1], f[2].upper(), f[3].upper(), f[4:]
    if len(gts) != n_samples:
        sys.exit(f"ERROR {len(gts)} genotypes at {pos}, expected {n_samples}")
    off = pos - 1 - win_start
    if not 0 <= off < input_len:
        sys.exit(f"ERROR variant at {pos} outside window")
    if ref_seq[off] != ref:
        sys.exit(f"ERROR REF mismatch at {pos}: VCF {ref} vs FASTA {ref_seq[off]}")
    v = len(variants)
    variants.append((pos, vid, ref, alt, off))
    for s, gt in enumerate(gts):
        phased = "|" in gt
        a = gt.replace("/", "|").split("|")
        if len(a) == 1:
            a = a * 2                            # haploid call -> same allele on both
        if "." in a:
            n_missing[s] += 1                    # missing allele -> left as REF
        al = [x == "1" for x in a]
        if al[0] != al[1] and not phased:
            n_unphased[s] += 1                   # het without phase: order taken as written
        for h in (0, 1):
            if al[h]:
                hap_vars[s][h].append(v)
        if pos == snp_pos:
            snp_dosage[s] = -1 if "." in a else int(al[0]) + int(al[1])
            snp_gt[s] = gt

print(f"SNVs in window: {len(variants)}; samples: {n_samples}")
print(f"Donors with unphased hets: {(n_unphased > 0).sum()}; with missing genotypes: {(n_missing > 0).sum()}")
if (snp_dosage < 0).all():
    sys.exit(f"ERROR focal SNP at {snp_pos} not found in genotype table")

# --- 6e. unique haplotypes: identical haplotypes are predicted once ---
hap_index = {(): 0}                              # key = tuple of ALT variant indices; () = REF
hap_ids = np.zeros((n_samples, 2), dtype=int)
for s in range(n_samples):
    for h in (0, 1):
        key = tuple(hap_vars[s][h])
        if key not in hap_index:
            hap_index[key] = len(hap_index)
        hap_ids[s, h] = hap_index[key]
keys = sorted(hap_index, key=hap_index.get)
print(f"Unique haplotypes to predict: {len(keys)} (including REF) for {2 * n_samples} donor haplotypes")

# --- 6f. one-hot encode: A=1000 C=0100 G=0010 T=0001 N=0000 ---
lut = {b: i for i, b in enumerate("ACGT")}
ref_oh = np.zeros((input_len, 4), dtype=np.float32)
for i, b in enumerate(ref_seq):
    if b in lut:
        ref_oh[i, lut[b]] = 1.0

def encode(key):
    x = ref_oh.copy()
    for v in key:
        off, alt = variants[v][4], variants[v][3]
        x[off] = 0.0
        x[off, lut[alt]] = 1.0
    return x

# --- 6g. predict unique haplotypes in batches ---
logcounts = np.zeros(len(keys))
coverage = np.zeros((len(keys), model.outputs[0].shape[1]), dtype=np.float32)
for b0 in range(0, len(keys), batch_size):
    batch = keys[b0:b0 + batch_size]
    x = np.stack([encode(k) for k in batch])
    profile, lc = model.predict(x, batch_size=len(batch), verbose=0)
    lc = np.asarray(lc).reshape(len(batch), -1)[:, 0]
    profile = np.asarray(profile)
    probs = np.exp(profile - profile.max(axis=1, keepdims=True))
    probs /= probs.sum(axis=1, keepdims=True)
    logcounts[b0:b0 + len(batch)] = lc
    coverage[b0:b0 + len(batch)] = probs * np.expm1(lc)[:, None]   # per-base predicted coverage
    print(f"  predicted {min(b0 + batch_size, len(keys))}/{len(keys)}")

# ChromBPNet's count head is trained on log(1 + total reads), so expm1 converts it back to reads
counts = np.expm1(logcounts)
ref_cnt = counts[0]

# --- 6h. write tables ---
meta = [peak, celltype, condition, fold]

with open(f"{outdir}/haplotype_predictions.tsv", "w") as fh:
    fh.write("peak_id\tcelltype\tcondition\tfold\tsample_id\thaplotype\thap_id\tn_alt\tlogcount_pred\tcount_pred\n")
    for s, sid in enumerate(samples):
        for h in (0, 1):
            k = hap_ids[s, h]
            fh.write("\t".join(meta + [sid, f"hap{h + 1}", str(k), str(len(keys[k])),
                                       f"{logcounts[k]:.6f}", f"{counts[k]:.4f}"]) + "\n")

dip = (counts[hap_ids[:, 0]] + counts[hap_ids[:, 1]]) / 2      # diploid = mean of the two haplotypes (counts)
with open(f"{outdir}/diploid_predictions.tsv", "w") as fh:
    fh.write("peak_id\tcelltype\tcondition\tfold\tsample_id\tpopulation\tsuperpopulation\t"
             "snp_gt\tsnp_dosage\thap1_counts\thap2_counts\tdiploid_counts\tref_counts\tn_unphased_het\tn_missing_gt\n")
    for s, sid in enumerate(samples):
        p, sp = pop.get(sid, ("NA", "NA"))
        fh.write("\t".join(meta + [sid, p, sp, snp_gt[s], str(snp_dosage[s]),
                                   f"{counts[hap_ids[s, 0]]:.4f}", f"{counts[hap_ids[s, 1]]:.4f}",
                                   f"{dip[s]:.4f}", f"{ref_cnt:.4f}",
                                   str(n_unphased[s]), str(n_missing[s])]) + "\n")

# variants in the window with ALT allele frequency across all donor haplotypes
alt_n = np.zeros(len(variants), dtype=int)
for s in range(n_samples):
    for h in (0, 1):
        for v in hap_vars[s][h]:
            alt_n[v] += 1
with open(f"{outdir}/variants_in_window.tsv", "w") as fh:
    fh.write("pos\tid\tref\talt\tdist_to_center\talt_haplotypes\talt_freq\n")
    for v, (pos, vid, ref, alt, off) in enumerate(variants):
        fh.write(f"{pos}\t{vid}\t{ref}\t{alt}\t{pos - 1 - center}\t{alt_n[v]}\t{alt_n[v] / (2 * n_samples):.5f}\n")

# per-base coverage of each unique haplotype; hap_ids maps donors to rows
np.savez_compressed(f"{outdir}/predicted_coverage.npz",
                    coverage=coverage, hap_ids=hap_ids, samples=np.array(samples))

# --- 6i. summary: predicted counts by focal SNP genotype ---
print(f"REF predicted counts: {ref_cnt:.1f}")
for d in (0, 1, 2):
    m = snp_dosage == d
    if m.any():
        print(f"SNP dosage {d}: n={m.sum()}  diploid_counts mean={dip[m].mean():.1f}  median={np.median(dip[m]):.1f}")
PY
    die "Prediction step failed"
fi

log "Done. Outputs in ${OUTDIR}/"


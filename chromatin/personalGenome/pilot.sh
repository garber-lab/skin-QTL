#!/usr/bin/env bash
#BSUB -J pilot_KRT_IFNgamma
#BSUB -q short
#BSUB -W 1:00
#BSUB -n 1
#BSUB -R "rusage[mem=8000]"
#BSUB -o pilot_KRT_IFNgamma.%J.log
#BSUB -e pilot_KRT_IFNgamma.%J.err

# pilot: haplotype-based personal chromBPNet prediction for 1 1000G donor and 1 ATAC-seq peak
#        using chromBPNet model KRT IFNgamma fold0
# steps:
# 1. SNP -> overlapping peak (strongest summit) -> 2114 bp model input window centered on the summit
# 2. Check the sample and chromosome exist in the VCF
# 3. Subset the phased 1000G VCF to the donor + window (biallelic SNVs, phased, non-missing)
# 4. Check VCF REF alleles against the FASTA
# 5. bcftools consensus -H 1 / -H 2 --> hap1 / hap2 sequences
# 6. Predict REF, hap1, hap2 with chrombpnet_nobias.h5
#
# Usage:
# bsub < pilot.sh

set -eo pipefail

BASE=/pi/manuel.garber-umw/human/skin/eQTLs
CHROMBPNET=${BASE}/chromatin/chrombpnet

# ---------------- settings ----------------
CELLTYPE=KRT
CONDITION=IFNgamma                        # filename-safe label written into output tables
MODEL_COND=IFN                            # label used in model/peak directory names
FOLD=0
SAMPLE=HG00107        # HG00107 is 0|1 at rs6844475 C>T
SNP=rs6844475         # strong eQTL in a KRT IFNG peak with strong chromBPNet scores
SNP_REF=C             # dbSNP REF allele, checked against hg38 below
INPUT_LEN=2114
RUN_MODEL=1           # 0 = stop after building sequences
PYTHON=python

MODEL=${CHROMBPNET}/ATAC_${CELLTYPE}_${MODEL_COND}_merged_alldonors/biasFactor0.5/fold${FOLD}/chrombpnet_model/models/chrombpnet_nobias.h5
PEAKS=${BASE}/ATACseq_DolphinNext/encode_atacseq_pipeline/MACS2_for_chrombpnet/ATAC_${CELLTYPE}_${MODEL_COND}_merged_alldonors_peaks.narrowPeak
FASTA=${CHROMBPNET}/ATAC_FRB_PBS_merged_alldonors/biasFactor0.5/data/downloads/hg38.fa
METADATA=${BASE}/literature/1000G/GRCh38_30x_phased/metadata/20130606_g1k_3202_samples_ped_population.txt
ALLSNPBED=${BASE}/DREG/2026models/data/all_genotyped_SNPs.bed
OUTDIR=${BASE}/chromatin/personalGenome/pilot/${SNP}_${SAMPLE}_${CELLTYPE}_${CONDITION}_fold${FOLD}
# -------------------------------------------

mkdir -p "${OUTDIR}"
LOG=${OUTDIR}/pilot.log
: > "${LOG}"
log() { echo "$(date '+%F %T') INFO $*" | tee -a "${LOG}"; }
die() { echo "$(date '+%F %T') ERROR $*" | tee -a "${LOG}" >&2; exit 1; }

module load samtools
module load bcftools
module load bedtools
set -u                                    # after module load: unset variables are now errors
cd "${OUTDIR}"

for f in "${PEAKS}" "${FASTA}" "${FASTA}.fai" "${ALLSNPBED}"; do
    [[ -e ${f} ]] || die "Missing file: ${f}"
done
if [[ ${RUN_MODEL} == 1 ]]; then [[ -e ${MODEL} ]] || die "Missing model: ${MODEL}"; fi

# ---------------- 1. SNP -> peak -> summit -> model input window ----------------

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



# ---------------- 2. sample checks in 1000VCF data, log metadata ----------------
bcftools query -l ${VCF} | grep -x ${SAMPLE} >/dev/null || die "Sample $SAMPLE not in VCF header"
bcftools view -h ${VCF} | grep "^##contig=<ID=$P_CHROM>" >/dev/null || die "$P_CHROM not in VCF contigs (chr1 vs 1?)"
POP=$(awk -v s=${SAMPLE} '$2 == s {print $6"\t"$7}' ${METADATA} 2>/dev/null || true)
log "Sample ${SAMPLE} population/superpopulation: ${POP:-NA}"

# ---------------- 3. variants in 2114bp window for donor of interest ----------------
# whole chr4 VCF, keep only sample and the 2114bp window -> get all variant record in the window with sample's genotype
ALL_VCF=${OUTDIR}/window_all.vcf.gz
SNV_VCF=${OUTDIR}/window_snvs.vcf.gz
bcftools view -s ${SAMPLE} -r ${REGION} -Oz -o ${ALL_VCF} ${VCF} 
bcftools index -f -t ${ALL_VCF}

# biallelic SNVs, non-missing, phased, sample carries >=1 ALT
bcftools view -v snps -m2 -M2 -p -i 'GT!="mis" && GT!="RR"' -Oz -o  ${SNV_VCF} ${ALL_VCF}
bcftools index -f -t ${SNV_VCF}

N_ALL=$(bcftools view -H ${ALL_VCF} | wc -l | tr -d ' ')
N_SNV=$(bcftools view -H ${SNV_VCF} | wc -l | tr -d ' ')

log "Variant records in window: $N_ALL; applied SNV sites (sample carries ALT): $N_SNV"
log "Genotype class counts (all records): $(bcftools query -f '[%GT]\n' "$ALL_VCF" | sort | uniq -c | awk '{printf "%s=%s ", $2, $1}')"
log "Skipped: indel=$(bcftools view -H -V snps "$ALL_VCF" | wc -l | tr -d ' ') \
multiallelic_snv=$(bcftools view -H -v snps -m3 "$ALL_VCF" | wc -l | tr -d ' ') \
missing_gt=$(bcftools view -H -i 'GT="mis"' "$ALL_VCF" | wc -l | tr -d ' ') \
unphased=$(bcftools view -H -P "$ALL_VCF" | wc -l | tr -d ' ')"

# check SNP
SNP_GT=$(bcftools query -r "${CHROM}:${SNP_POS}" -f '[%GT]' "${SNV_VCF}")
[[ -n ${SNP_GT} ]] || die "${SNP} (${CHROM}:${SNP_POS}) not among applied SNVs for ${SAMPLE}"
log "${SNP} genotype in ${SAMPLE}: ${SNP_GT}"

# ---------------- 4. REF check ----------------
# bcftools norm --check-ref e exits with an error on any REF/FASTA mismatch
bcftools norm -f "$FASTA" --check-ref e -Ou "$ALL_VCF" > /dev/null 2>"$OUTDIR/refcheck.err" \
    || die "REF allele mismatch vs FASTA (see $OUTDIR/refcheck.err)"
log "REF alleles match FASTA"

bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT]\n" "$SNV_VCF" \
    | awk -v c="$CENTER" 'BEGIN{OFS="\t"; print "chrom","pos","id","ref","alt","gt","dist_to_center"} {print $0, $2 - 1 - c}' \
    > "$OUTDIR/variants_applied.tsv"

# ---------------- 5. haplotype sequences ----------------
# build three 2114bp sequences (REF, hap1, hap2), saves to one FASTA file

# 5a. helper function: turns a FASTA file into one plain uppercase string
seq_of() { grep -v '^>' "$1" | tr -d '\n' | tr 'acgtn' 'ACGTN'; }

# 5b. builds 3 sequences.
# bcftools consensus writes the donor's alleles into the window.
samtools faidx ${FASTA} ${REGION} > ${OUTDIR}/REF.fa
for H in 1 2; do
    samtools faidx "$FASTA" "$REGION" | bcftools consensus -s "$SAMPLE" -H "$H" "$SNV_VCF" > "$OUTDIR/hap$H.fa" 2>>"$LOG"
done

# 5c. Put 3 seq into one FASTA and check
REF_SEQ=$(seq_of "$OUTDIR/REF.fa")
SEQ_FA=$OUTDIR/sequences.fa
: > "$SEQ_FA"  # start with empty file
for LABEL in REF hap1 hap2; do
    S=$(seq_of "$OUTDIR/$LABEL.fa")
    N_DIFF=$(paste <(echo "$S" | fold -w1) <(echo "$REF_SEQ" | fold -w1) | awk '$1 != $2' | wc -l | tr -d ' ')
    log "$LABEL: length ${#S}, $N_DIFF bases differ from reference"
    printf ">%s|%s|%s|%s\n%s\n" "$P_NAME" "$SAMPLE" "$LABEL" "$REGION" "$S" >> "$SEQ_FA"  	# append header + sequence
done
rm -f "$OUTDIR/REF.fa" "$OUTDIR/hap1.fa" "$OUTDIR/hap2.fa"


# 5d. Check: ALT alleles in = bases changed out
for H in 1 2; do
    # EXP = n of ALT alleles hap H should carry
    EXP=$(awk -F'\t' -v h="$H" 'NR > 1 { split($6, g, "|"); if (g[h] != "0") n++ } END {print n + 0}' "$OUTDIR/variants_applied.tsv")
    OBS=$(awk -v h="$H" 'NR == 2*(h+1) { print }' "$SEQ_FA" | fold -w1 | paste - <(echo "$REF_SEQ" | fold -w1) | awk '$1 != $2' | wc -l | tr -d ' ')
    [[ $EXP == "$OBS" ]] || die "hap$H: expected $EXP ALT bases but sequence differs at $OBS positions"
done
log "Haplotype QC passed (ALT counts match, 0|0 reproduces reference)"




# ---------------- 6. ChromBPNet prediction (inside the chrombpnet container) ----------------
# The Python code between <<'PY' and PY is fed straight into the container's python and run here.
# Arguments after "python -" arrive in Python as sys.argv[1:].
SIF="/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif"

log "Predicting with ${MODEL}"
if ! singularity exec "${SIF}" python - \
        "${MODEL}" "${SEQ_FA}" "${OUTDIR}" "${INPUT_LEN}" \
        "${P_NAME}" "${SAMPLE}" "${CELLTYPE}" "${CONDITION}" "${FOLD}" \
        <<'PY' 2>&1 | tee -a "${LOG}"; then
import sys
import numpy as np
import tensorflow as tf

# --- 6a. arguments from bash ---
model_path, seq_fa, outdir, input_len, peak, sample, celltype, condition, fold = sys.argv[1:]
input_len = int(input_len)
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

# --- 6c. read sequences.fa: header label (REF/hap1/hap2) + sequence ---
labels, seqs = [], []
for line in open(seq_fa):
    if line.startswith(">"):
        labels.append(line[1:].split("|")[2])
    else:
        seqs.append(line.strip())

# --- 6d. one-hot encode: A=1000 C=0100 G=0010 T=0001 N=0000 ---
x = np.zeros((len(seqs), input_len, 4), dtype=np.float32)
for k, s in enumerate(seqs):
    for i, b in enumerate(s):
        j = "ACGT".find(b)
        if j >= 0:
            x[k, i, j] = 1.0

# --- 6e. predict: profile logits (shape over central 1000 bp) + logcounts (total signal) ---
profile, logcounts = model.predict(x, batch_size=len(seqs), verbose=0)
logcounts = np.asarray(logcounts).reshape(len(seqs), -1)[:, 0]
profile = np.asarray(profile)

# ChromBPNet's count head is trained on log(1 + total reads), so expm1 converts it back to reads
counts = np.expm1(logcounts)

# per-base predicted coverage over the central 1000 bp = softmax(profile logits) * total counts
probs = np.exp(profile - profile.max(axis=1, keepdims=True))
probs /= probs.sum(axis=1, keepdims=True)
coverage = probs * counts[:, None]

pred_log = dict(zip(labels, map(float, logcounts)))
pred_cnt = dict(zip(labels, map(float, counts)))
cov = dict(zip(labels, coverage))

# diploid = average of the two haplotypes, on the counts scale (not log)
dip_cnt = (pred_cnt["hap1"] + pred_cnt["hap2"]) / 2
dip_cov = (cov["hap1"] + cov["hap2"]) / 2

np.savez_compressed(f"{outdir}/predicted_coverage.npz",
                    REF=cov["REF"], hap1=cov["hap1"], hap2=cov["hap2"], diploid=dip_cov)

# --- 6f. write tables ---
meta = [peak, sample, celltype, condition, fold]
with open(f"{outdir}/haplotype_predictions.tsv", "w") as fh:
    fh.write("peak_id\tsample_id\tcelltype\tcondition\tfold\thaplotype\tlogcount_pred\tcount_pred\n")
    for lab in labels:
        fh.write("\t".join(meta + [lab, f"{pred_log[lab]:.6f}", f"{pred_cnt[lab]:.4f}"]) + "\n")

with open(f"{outdir}/diploid_predictions.tsv", "w") as fh:
    fh.write("peak_id\tsample_id\tcelltype\tcondition\tfold\thap1_counts\thap2_counts\tdiploid_counts\tref_counts\n")
    fh.write("\t".join(meta + [f"{pred_cnt['hap1']:.4f}", f"{pred_cnt['hap2']:.4f}",
                               f"{dip_cnt:.4f}", f"{pred_cnt['REF']:.4f}"]) + "\n")

print(f"predicted counts REF={pred_cnt['REF']:.1f} hap1={pred_cnt['hap1']:.1f} hap2={pred_cnt['hap2']:.1f} "
      f"diploid={dip_cnt:.1f} hap2/hap1={pred_cnt['hap2'] / pred_cnt['hap1']:.3f}")

PY
    die "Prediction step failed"
fi

log "Done. Outputs in ${OUTDIR}/"

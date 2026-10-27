#!/bin/bash
#BSUB -J fill_A2_ssc
#BSUB -q short
#BSUB -W 2:00
#BSUB -n 1
#BSUB -R "rusage[mem=8000]"
#BSUB -o /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/log/fill_A2_%J.out
#BSUB -e /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/log/fill_A2_%J.err
# ---------------------------------------------------------------------------
# Fill in A2 for the Lopez-Isac SSc meta-analysis using the hg19 reference.
#
# WHAT A REFERENCE GENOME CAN AND CANNOT TELL YOU
#
#   hg19.fa gives the REFERENCE base at a position. It does not know the
#   variant's allele pair. So:
#
#     A1 != REF  ->  the file's A1 is the alternate allele, therefore
#                    A2 = REF.  RESOLVED.
#     A1 == REF  ->  A1 is the reference allele, so A2 is the ALT -- and the
#                    reference genome has no idea what the ALT is.
#                    NOT RESOLVABLE HERE. Left as '?' and listed in the
#                    .unresolved file for a dbSNP/1000G lookup by rsID.
#
#   PLINK sets A1 to the MINOR allele by default, and for most variants the
#   minor allele is the non-reference one, so the resolvable fraction should be
#   large. The QC block prints the actual split -- do not assume it.
#
# Input is read in file order and the reference is loaded one chromosome at a
# time, so memory stays at ~one chromosome (chr1, ~250 MB) regardless of how
# many variants there are. This requires each chromosome's rows to be
# contiguous, which they are; the script aborts if that is violated rather
# than silently thrashing.
#
# Submit:  bsub < fill_A2_from_hg19.sh
# ---------------------------------------------------------------------------
set -euo pipefail

dir=${COLOC_DIR:-/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc}
trait=${TRAIT:-systemic_sclerosis_lopez2019}
IN=${IN:-Lopez-Isac_META_GWAS_SSc.noN.tsv}
OUT=${OUT:-Lopez-Isac_META_GWAS_SSc.A2filled.tsv}
REF=${REF:-/share/data/umw_biocore/genome_data/human/hg19/hg19.fa}

cd "${dir}/${trait}"
[ -s "$IN" ]        || { echo "missing input: $(pwd)/$IN" >&2; exit 1; }
[ -s "$REF" ]       || { echo "missing reference: $REF" >&2; exit 1; }
[ -s "${REF}.fai" ] || { echo "missing index: ${REF}.fai  (samtools faidx $REF)" >&2; exit 1; }

python3 - "$IN" "$OUT" "$REF" <<'PY'
import sys, os
infile, outfile, reffa = sys.argv[1], sys.argv[2], sys.argv[3]
faifile = reffa + ".fai"

# ---- .fai: name -> (length, byte offset, bases per line, bytes per line) ----
fai = {}
for line in open(faifile):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 5:
        fai[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))

def fai_key(c):
    """The meta file uses 1..22/X/Y; UCSC hg19.fa uses chr1..chrY. Try both."""
    for k in ("chr" + c, c, "chr" + c.upper(), c.upper()):
        if k in fai:
            return k
    return None

def load_chrom(key):
    """Read one chromosome into memory, uppercased (hg19.fa is soft-masked)."""
    length, offset, linebases, linewidth = fai[key]
    nbytes = (length // linebases) * linewidth + (length % linebases)
    with open(reffa, "rb") as fh:
        fh.seek(offset)
        raw = fh.read(nbytes)
    return raw.replace(b"\n", b"").replace(b"\r", b"").upper()

fin  = open(infile)
fout = open(outfile, "w")
funr = open(outfile + ".unresolved", "w")

hdr = fin.readline().rstrip("\n").split()
ix  = {name: i for i, name in enumerate(hdr)}
for need in ("CHR", "BP", "SNP", "A1", "A2"):
    if need not in ix:
        sys.exit(f"ERROR: no '{need}' column in {infile}\nheader: {hdr}")
C, B, S, A1, A2 = ix["CHR"], ix["BP"], ix["SNP"], ix["A1"], ix["A2"]

fout.write("\t".join(hdr) + "\n")
funr.write("CHR\tBP\tSNP\tA1\treason\n")

seq = None; cur = None; done = set()
n = resolved = same_as_ref = indel = nocontig = refN = 0
a1len = {}

for line in fin:
    f = line.split()
    if len(f) < len(hdr):
        continue
    n += 1
    c, bp, a1 = f[C], int(f[B]), f[A1].upper()
    a1len[len(a1)] = a1len.get(len(a1), 0) + 1

    if c != cur:
        if c in done:
            sys.exit(f"ERROR: chromosome {c} reappears after another chromosome.\n"
                     f"       This script streams the reference one chromosome at a time and\n"
                     f"       needs each chromosome's rows contiguous. Sort the input first:\n"
                     f"       (head -1 {infile}; tail -n +2 {infile} | sort -k1,1V -k2,2n) > sorted.tsv")
        if cur is not None:
            done.add(cur)
        cur = c
        key = fai_key(c)
        seq = load_chrom(key) if key else None
        print(f"  chr{c}: {'loaded ' + key if key else 'NOT IN REFERENCE - skipping'}",
              file=sys.stderr, flush=True)

    if seq is None:
        nocontig += 1
        funr.write(f"{c}\t{bp}\t{f[S]}\t{a1}\tcontig_not_in_reference\n")
        fout.write("\t".join(f[:len(hdr)]) + "\n")
        continue

    if len(a1) != 1:
        indel += 1
        funr.write(f"{c}\t{bp}\t{f[S]}\t{a1}\tnot_a_single_base_allele\n")
        fout.write("\t".join(f[:len(hdr)]) + "\n")
        continue

    ref = chr(seq[bp - 1]) if 0 < bp <= len(seq) else "N"
    if ref == "N":
        refN += 1
        funr.write(f"{c}\t{bp}\t{f[S]}\t{a1}\treference_base_is_N\n")
    elif ref != a1:
        f[A2] = ref                      # A1 is the ALT, so the other allele is REF
        resolved += 1
    else:
        same_as_ref += 1                 # A1 IS the reference; ALT unknown from a FASTA
        funr.write(f"{c}\t{bp}\t{f[S]}\t{a1}\tA1_equals_reference_need_dbsnp\n")

    fout.write("\t".join(f[:len(hdr)]) + "\n")

for fh in (fin, fout, funr):
    fh.close()

pct = lambda x: 100.0 * x / n if n else 0.0
print(f"""
================= A2 fill summary =================
  variants processed          {n:>12,}
  A2 RESOLVED (A1 != REF)     {resolved:>12,}  ({pct(resolved):5.2f}%)
  A2 still '?' :
    A1 == REF (need dbSNP)    {same_as_ref:>12,}  ({pct(same_as_ref):5.2f}%)
    non-single-base A1        {indel:>12,}  ({pct(indel):5.2f}%)
    reference base is N       {refN:>12,}  ({pct(refN):5.2f}%)
    contig not in reference   {nocontig:>12,}  ({pct(nocontig):5.2f}%)

  A1 allele lengths: {dict(sorted(a1len.items()))}

  BUILD SANITY CHECK -- read 'A1 == REF' above:
    PLINK's A1 is the minor allele, which is usually the non-reference one, so
    on correct hg19 coordinates this should land roughly in the 5-25% range.
    A value near 25% with no clear peak is what random positions look like:
    if you see that, suspect the coordinates are not hg19 after all.
    A value near 0% or near 100% means A1/A2 are not what this script assumes.

  unresolved list -> {os.path.basename(outfile)}.unresolved
===================================================""", file=sys.stderr)
PY

echo "done: $(pwd)/${OUT}"
head -3 "$OUT"
echo "unresolved (first 3):"
head -4 "${OUT}.unresolved"

# ---------------------------------------------------------------------------
# NEXT STEP for the rows still marked '?': they need a variant database keyed by
# rsID, not a reference genome. dbSNP or the 1000G Phase III panel the study was
# imputed against will give the full allele pair. Do that lookup BEFORE the
# hg19 -> hg38 liftOver, while these coordinates still match the panel.
# ---------------------------------------------------------------------------

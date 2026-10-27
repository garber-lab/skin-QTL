#!/usr/bin/env python3
"""
Fill in A2 for the Lopez-Isac 2019 SSc meta-analysis from the 1000 Genomes
Phase 3 sites VCF (GRCh37).

WHY A PANEL AND NOT A REFERENCE GENOME
    hg19.fa stores one base per position -- the REFERENCE allele. It can only
    resolve A2 when A1 != REF (then A2 = REF). When A1 IS the reference base,
    the alternate allele is simply not in a FASTA, and those variants stay '?'.

    The 1000G Phase 3 panel stores the allele PAIR, so it resolves both cases.
    It is also the correct panel by construction: the paper imputed against it
    (Lopez-Isac 2019 Nat Commun, PMC6823490, Methods -- "imputation using
    IMPUTE2 and the 1000 Genome Project Phase III data as reference panel"),
    so every variant in the file came out of this panel.

    Two further things it gives that a FASTA cannot:
      - EUR_AF, which this study deposited nowhere. Finemapping needs a MAF.
      - a check that A1 is actually one of the site's two alleles, which
        catches strand problems and coordinate drift. A FASTA cannot tell
        "A1 is the ALT" from "A1 is wrong".

MEMORY
    One dict keyed by packed (chrom, pos) over the ~4.7M meta variants, filled
    in a single streaming pass over the VCF. ~1-2 GB. The VCF is never held in
    memory and the meta file is re-read from disk for the output pass.

Usage:
    fill_A2_from_1000G.py <meta.txt> <1000G.sites.vcf.gz> <out.tsv>
"""
import sys, gzip, os

if len(sys.argv) != 4:
    sys.exit(__doc__)
metafile, vcffile, outfile = sys.argv[1:4]

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
revcomp = lambda s: "".join(COMP.get(b, "N") for b in reversed(s))

def chrom_idx(c):
    c = c.lstrip("chr").lstrip("CHR")
    if c == "X" or c == "23": return 23
    if c == "Y" or c == "24": return 24
    if c in ("MT", "M"):      return 25
    return int(c) if c.isdigit() else None

def key_of(c, pos):
    i = chrom_idx(c)
    return None if i is None else (i << 32) | pos

# ---------------------------------------------------------------- pass 1 ----
# Which positions do we actually need? Everything else in the VCF is skipped.
sys.stderr.write("pass 1/3: reading meta positions ...\n")
with open(metafile) as fh:
    hdr = fh.readline().split()
    ix = {n: i for i, n in enumerate(hdr)}
    for need in ("CHR", "BP", "SNP", "A1", "A2"):
        if need not in ix:
            sys.exit(f"ERROR: no '{need}' column. header={hdr}")
    C, B, S, A1C, A2C = ix["CHR"], ix["BP"], ix["SNP"], ix["A1"], ix["A2"]
    want = {}
    for line in fh:
        f = line.split()
        if len(f) < len(hdr):
            continue
        k = key_of(f[C], int(f[B]))
        if k is not None:
            want[k] = None
sys.stderr.write(f"  {len(want):,} distinct positions wanted\n")

# ---------------------------------------------------------------- pass 2 ----
# Stream the sites VCF once. Parse the cheap prefix first and bail out before
# touching INFO for the ~95% of lines we do not care about.
sys.stderr.write("pass 2/3: streaming 1000G sites VCF ...\n")
seen = nline = 0
with gzip.open(vcffile, "rt") as fh:
    for line in fh:
        if line[0] == "#":
            continue
        nline += 1
        if nline % 10_000_000 == 0:
            sys.stderr.write(f"  {nline:,} VCF lines, {seen:,} matched\n"); sys.stderr.flush()
        t1 = line.find("\t")
        t2 = line.find("\t", t1 + 1)
        k = key_of(line[:t1], int(line[t1 + 1:t2]))
        if k is None or k not in want:
            continue
        f = line.rstrip("\n").split("\t")
        ref, alts = f[3].upper(), f[4].upper().split(",")
        eur = []
        for kv in f[7].split(";"):
            if kv.startswith("EUR_AF="):
                eur = kv[7:].split(","); break
        rec = want[k]
        if rec is None:
            rec = want[k] = []
        for i, alt in enumerate(alts):
            rec.append((ref, alt, eur[i] if i < len(eur) else ""))
        seen += 1
sys.stderr.write(f"  {nline:,} VCF lines scanned, {seen:,} records matched\n")

# ---------------------------------------------------------------- pass 3 ----
sys.stderr.write("pass 3/3: resolving A2 and writing ...\n")
fin  = open(metafile)
fout = open(outfile, "w")
funr = open(outfile + ".unresolved", "w")
fin.readline()

out_hdr = list(hdr) + ["EUR_AF_1000G"]
fout.write("\t".join(out_hdr) + "\n")
funr.write("CHR\tBP\tSNP\tA1\treason\n")

n = res_alt = res_ref = 0
notfound = mismatch = ambig = strand = 0

for line in fin:
    f = line.split()
    if len(f) < len(hdr):
        continue
    n += 1
    c, bp, snp, a1 = f[C], int(f[B]), f[S], f[A1C].upper()
    recs = want.get(key_of(c, bp))
    af = ""

    if not recs:
        notfound += 1
        funr.write(f"{c}\t{bp}\t{snp}\t{a1}\tposition_not_in_1000G\n")
    else:
        # A1 may be the ALT (A2=REF) or the REF (A2=ALT). Collect every
        # consistent answer; a unique one is a resolution, several is ambiguity.
        as_alt = {(r, a, e) for (r, a, e) in recs if a1 == a}
        as_ref = {(r, a, e) for (r, a, e) in recs if a1 == r}
        cand = {r for (r, a, e) in as_alt} | {a for (r, a, e) in as_ref}
        if len(cand) == 1:
            f[A2C] = cand.pop()
            if as_alt:
                res_alt += 1
                e = [e for (_, _, e) in as_alt if e]
                af = e[0] if e else ""
            else:
                res_ref += 1
                tot = sum(float(e) for (_, _, e) in as_ref if e)
                af = f"{1.0 - tot:.6g}" if any(e for (_, _, e) in as_ref) else ""
        elif len(cand) > 1:
            ambig += 1
            funr.write(f"{c}\t{bp}\t{snp}\t{a1}\tmultiple_possible_A2:{'/'.join(sorted(cand))}\n")
        else:
            rc = revcomp(a1)
            if any(rc == a or rc == r for (r, a, e) in recs):
                strand += 1
                funr.write(f"{c}\t{bp}\t{snp}\t{a1}\tA1_matches_only_on_opposite_strand\n")
            else:
                mismatch += 1
                obs = "/".join(sorted({f"{r}>{a}" for (r, a, e) in recs}))
                funr.write(f"{c}\t{bp}\t{snp}\t{a1}\tA1_not_an_allele_here:{obs}\n")

    fout.write("\t".join(f[:len(hdr)] + [af if af else "NA"]) + "\n")

for fh in (fin, fout, funr):
    fh.close()

pct = lambda x: 100.0 * x / n if n else 0.0
resolved = res_alt + res_ref
sys.stderr.write(f"""
=================== A2 fill summary ===================
  variants processed            {n:>12,}
  A2 RESOLVED                   {resolved:>12,}  ({pct(resolved):6.2f}%)
      A1 was the ALT (A2=REF)   {res_alt:>12,}  ({pct(res_alt):6.2f}%)
      A1 was the REF (A2=ALT)   {res_ref:>12,}  ({pct(res_ref):6.2f}%)   <-- a FASTA CANNOT do these
  still '?' :
      position not in 1000G     {notfound:>12,}  ({pct(notfound):6.2f}%)
      A1 not an allele there    {mismatch:>12,}  ({pct(mismatch):6.2f}%)
      ambiguous (multi-allelic) {ambig:>12,}  ({pct(ambig):6.2f}%)
      matches only on rev-comp  {strand:>12,}  ({pct(strand):6.2f}%)

  HOW TO READ THIS
    'A1 not an allele there' should be near zero. This file was imputed FROM
    this panel, so a real allele mismatch means the coordinates or the alleles
    are not what we think. A large 'matches only on rev-comp' would mean the
    deposit is on the opposite strand from the panel.

  unresolved list -> {os.path.basename(outfile)}.unresolved
=======================================================\n""")

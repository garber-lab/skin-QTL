#!/usr/bin/env python3
"""
Lopez-Isac 2019 SSc: drop unresolved variants, derive beta/se from OR and P,
and lift hg19 -> hg38.

ORDER MATTERS. The liftOver is last on purpose: it is the step that invalidates
the GRCh37 position join used to fill A2 and EUR_AF, so everything that needs
GRCh37 coordinates has to happen first.

BETA / SE
    beta = ln(OR)
    z    = -Phi^-1(P/2)          two-sided, via statistics.NormalDist (stdlib,
                                 Wichura AS241; exact at P=0.05 -> 1.959963985,
                                 and fine at this file's min P of 3.25e-43)
    se   = |beta| / z
    varbeta = se^2

    Rows that cannot yield a usable se are DROPPED, not fudged:
      OR == 1  -> beta == 0 -> se == 0 -> varbeta == 0, which coloc cannot use
      P  == 1  -> z == 0    -> se undefined (0/0)
    Both are real in this file (4,527 and 187 rows).

    PRECISION CAVEAT: se is only as precise as the deposited P, which carries
    ~4 significant figures. Round-tripping a known se through a 4-sig-fig P and
    back recovers it to ~1e-3 relative at very small P, exactly at moderate P.
    That is inherent to a deposit with no SE column, not a defect here.

LIFTOVER
    pyliftover + UCSC hg19ToHg38.over.chain.gz. pyliftover is 0-based, the file
    is 1-based, so positions are shifted on the way in and back on the way out.

    Dropped: unmapped positions, positions with more than one target (genuinely
    ambiguous), and positions that land on a different chromosome.

    STRAND: a chain block can be on the minus strand, and UCSC liftOver does not
    touch alleles. Where that happens A1/A2 are reverse-complemented, otherwise
    the alleles would silently disagree with the new coordinates. The count is
    reported.

    The result is re-sorted: hg19 order does not imply hg38 order across
    rearranged regions, and anything downstream wanting tabix needs sorted input.

Usage:
    ssc_beta_se_liftover.py <in.tsv.gz> <chain.gz> <out.tsv>
"""
import sys, gzip, math, os
from statistics import NormalDist

if len(sys.argv) != 4:
    sys.exit(__doc__)
infile, chainfile, outfile = sys.argv[1:4]

from pyliftover import LiftOver
nd = NormalDist()
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
revcomp = lambda s: "".join(COMP.get(b, "N") for b in reversed(s))

sys.stderr.write("loading chain ...\n")
lo = LiftOver(chainfile)

opener = gzip.open if infile.endswith(".gz") else open
fin = opener(infile, "rt")
hdr = fin.readline().rstrip("\n").split("\t")
ix = {n: i for i, n in enumerate(hdr)}
for need in ("CHR", "BP", "SNP", "A1", "A2", "P", "OR"):
    if need not in ix:
        sys.exit(f"ERROR: no '{need}' column. header={hdr}")
C, B, S, A1, A2, P, OR = (ix[k] for k in ("CHR", "BP", "SNP", "A1", "A2", "P", "OR"))
AF = ix.get("EUR_AF_1000G"); Q = ix.get("Q"); I = ix.get("I")

OUT = ["CHR", "BP", "SNP", "A1", "A2", "BETA", "SE", "VARBETA", "P", "EUR_AF", "Q", "I"]
rows = []
n = 0
d_ambig = d_or = d_p = d_beta0 = d_z0 = d_se = 0
d_unmapped = d_multi = d_xchrom = 0
flipped = 0

sys.stderr.write("converting ...\n")
for line in fin:
    f = line.rstrip("\n").split("\t")
    if len(f) < len(hdr):
        continue
    n += 1
    if n % 1_000_000 == 0:
        sys.stderr.write(f"  {n:,} rows\n"); sys.stderr.flush()

    a1, a2 = f[A1].upper(), f[A2].upper()
    if a2 == "?" or a1 == "?":
        d_ambig += 1; continue

    # ---- beta / se -------------------------------------------------------
    try:
        orv = float(f[OR]); pv = float(f[P])
    except ValueError:
        d_or += 1; continue
    if not (orv > 0) or not math.isfinite(orv):
        d_or += 1; continue
    if not (0 < pv <= 1):
        d_p += 1; continue
    beta = math.log(orv)
    if beta == 0.0:                      # OR exactly 1 -> se would be 0
        d_beta0 += 1; continue
    z = -nd.inv_cdf(pv / 2.0)
    if z <= 0:                           # P exactly 1
        d_z0 += 1; continue
    se = abs(beta) / z
    if not math.isfinite(se) or se <= 0:
        d_se += 1; continue

    # ---- liftOver --------------------------------------------------------
    src = "chr" + f[C]
    res = lo.convert_coordinate(src, int(f[B]) - 1)
    if not res:
        d_unmapped += 1; continue
    if len(res) > 1:
        d_multi += 1; continue
    tchrom, tpos0, tstrand = res[0][0], res[0][1], res[0][2]
    if tchrom != src:
        d_xchrom += 1; continue
    if tstrand == "-":
        a1, a2 = revcomp(a1), revcomp(a2)
        flipped += 1

    rows.append((tchrom[3:], tpos0 + 1, f[S], a1, a2,
                 # .9g not .6g: at 6 sig figs the printed varbeta and the square
                 # of the printed se disagree by up to 1.5e-5 relative, which is
                 # harmless for coloc but looks like an inconsistency on audit.
                 f"{beta:.9g}", f"{se:.9g}", f"{se*se:.9g}", f[P],
                 f[AF] if AF is not None else "NA",
                 f[Q] if Q is not None else "NA",
                 f[I] if I is not None else "NA"))
fin.close()

# hg19 order does not imply hg38 order; sort so the file is tabix-ready.
sys.stderr.write(f"sorting {len(rows):,} rows by hg38 coordinate ...\n")
def ckey(c):
    return int(c) if c.isdigit() else {"X": 23, "Y": 24, "MT": 25, "M": 25}.get(c, 99)
rows.sort(key=lambda r: (ckey(r[0]), r[1]))

with open(outfile, "w") as fh:
    fh.write("\t".join(OUT) + "\n")
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")

kept = len(rows)
pct = lambda x: 100.0 * x / n if n else 0.0
sys.stderr.write(f"""
============ beta/se + liftOver summary ============
  input rows                    {n:>12,}
  OUTPUT ROWS                   {kept:>12,}  ({pct(kept):6.2f}%)

  dropped -- unresolved A2      {d_ambig:>12,}  ({pct(d_ambig):6.2f}%)
  dropped -- beta/se not usable:
      OR == 1 (varbeta would=0) {d_beta0:>12,}  ({pct(d_beta0):6.2f}%)
      P  == 1 (z == 0)          {d_z0:>12,}  ({pct(d_z0):6.2f}%)
      OR unparseable/<=0        {d_or:>12,}  ({pct(d_or):6.2f}%)
      P out of (0,1]            {d_p:>12,}  ({pct(d_p):6.2f}%)
      se non-finite             {d_se:>12,}  ({pct(d_se):6.2f}%)
  dropped -- liftOver:
      unmapped in hg38          {d_unmapped:>12,}  ({pct(d_unmapped):6.2f}%)
      maps to >1 location       {d_multi:>12,}  ({pct(d_multi):6.2f}%)
      maps to other chromosome  {d_xchrom:>12,}  ({pct(d_xchrom):6.2f}%)

  alleles reverse-complemented (minus-strand chain block): {flipped:,}
====================================================\n""")

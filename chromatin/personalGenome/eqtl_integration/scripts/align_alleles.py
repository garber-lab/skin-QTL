"""
eQTL integration: align the eQTL genotypes (old VCF, REF not fixed to hg38) to the phased hg38 VCF used by ChromBPNet,
per chromosome. Run by align_alleles.sh (inside the container).

  python align_alleles.py <old.tsv.gz> <old_samples.txt> <new.tsv.gz> <new_samples.txt> <out.tsv.gz>

Inputs are `bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n'` and `bcftools query -l` of each VCF.
old = skineQTL_imputed_hg38_filtered.vcf.gz (what the eQTL used; beta = per copy of its ALT allele)
new = phased/chr<N>_imputed_combined_hg38_rsID_filtered_phased.vcf.gz (REF = hg38; ChromBPNet ALT allele)

Donor IDs are reduced to F<number> the same way as the eQTL pipeline (step0.3): skineQTL-F064_skineQTL-F064 -> F64,
F25_F25 -> F25.

Per SNP (matched by rsID and position), two independent checks:
  allele_match  same | swap | strand | strand_swap | mismatch   (REF/ALT strings; A/T and C/G SNPs are ambiguous)
  r_dosage      correlation of old-ALT vs new-ALT dosage across the shared donors (+1 same, -1 swapped)
  flip          +1: old ALT = new ALT -> beta keeps its sign; -1: swapped -> multiply beta by -1; NaN: unresolved
                decided by r_dosage when |r| >= 0.8 (and n_shared >= 10), else by the alleles (non-palindromic only).
  agree         allele-based and genotype-based orientation agree (False = check by hand)
Only the beta orientation is needed downstream: aligned beta = beta * flip = effect per copy of the hg38 ALT allele,
the same allele as the ChromBPNet prediction.
"""
import sys

import numpy as np
import pandas as pd

COMP = str.maketrans("ACGT", "TGCA")
GT = {"0/0": 0, "0|0": 0, "0/1": 1, "1/0": 1, "0|1": 1, "1|0": 1, "1/1": 2, "1|1": 2}


def norm_id(s):
    """skineQTL-F064_skineQTL-F064 -> F64; F25_F25 -> F25 (same as the eQTL pipeline)."""
    s = s.split("_")[0].replace("skineQTL-", "")
    return "F" + str(int(s[1:])) if s[:1] == "F" and s[1:].isdigit() else s


def load(tsv, samples):
    ids = [norm_id(x.strip()) for x in open(samples) if x.strip()]
    d = pd.read_csv(tsv, sep="\t", header=None, dtype=str,
                    names=["chrom", "pos", "rsid", "ref", "alt"] + ids)
    d = d[~d.alt.str.contains(",")]                                   # biallelic only
    dos = d[ids].apply(lambda c: c.map(GT)).to_numpy(float)          # ALT dosage, NaN = missing / other
    return d[["chrom", "pos", "rsid", "ref", "alt"]].reset_index(drop=True), dos, ids


def allele_match(oR, oA, nR, nA):
    if (oR, oA) == (nR, nA):
        return "same"
    if (oR, oA) == (nA, nR):
        return "swap"
    if (oR.translate(COMP), oA.translate(COMP)) == (nR, nA):
        return "strand"
    if (oR.translate(COMP), oA.translate(COMP)) == (nA, nR):
        return "strand_swap"
    return "mismatch"


def rowcorr(a, b):
    """Pearson r per row over columns where both are present."""
    ok = ~np.isnan(a) & ~np.isnan(b)
    n = ok.sum(1)
    a, b = np.where(ok, a, 0), np.where(ok, b, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        ma, mb = a.sum(1) / n, b.sum(1) / n
        cov = (a * b).sum(1) / n - ma * mb
        va, vb = (a * a).sum(1) / n - ma ** 2, (b * b).sum(1) / n - mb ** 2
        return cov / np.sqrt(va * vb), n


def main(old_tsv, old_s, new_tsv, new_s, out):
    o, odos, oid = load(old_tsv, old_s)
    n, ndos, nid = load(new_tsv, new_s)
    shared = [x for x in nid if x in oid]
    print(f"donors: old {len(oid)}, new {len(nid)}, shared {len(shared)}")
    print(f"  only in old (eQTL, not phased): {sorted(set(oid) - set(nid))}")
    print(f"  only in new: {sorted(set(nid) - set(oid))}")

    o["oi"], n["ni"] = np.arange(len(o)), np.arange(len(n))
    m = o.merge(n, on=["chrom", "pos", "rsid"], suffixes=("_old", "_new"))
    oc, nc = [oid.index(x) for x in shared], [nid.index(x) for x in shared]
    r, ns = rowcorr(odos[m.oi.values][:, oc], ndos[m.ni.values][:, nc])

    m["allele_match"] = [allele_match(*x) for x in zip(m.ref_old, m.alt_old, m.ref_new, m.alt_new)]
    m["palindromic"] = [{a, b} in ({"A", "T"}, {"C", "G"}) for a, b in zip(m.ref_new, m.alt_new)]
    m["r_dosage"], m["n_shared"] = r, ns
    by_allele = m.allele_match.map({"same": 1, "strand": 1, "swap": -1, "strand_swap": -1}).astype(float)
    by_allele[m.palindromic] = np.nan
    by_geno = pd.Series(np.where((np.abs(r) >= 0.8) & (ns >= 10), np.sign(r), np.nan), index=m.index)
    m["flip"] = by_geno.fillna(by_allele)
    m["agree"] = (by_allele.isna() | by_geno.isna()) | (by_allele == by_geno)
    m["old_alt_freq"] = np.nanmean(odos[m.oi.values], axis=1) / 2
    m["new_alt_freq"] = np.nanmean(ndos[m.ni.values], axis=1) / 2
    m = m.drop(columns=["oi", "ni"])
    m.to_csv(out, sep="\t", index=False, float_format="%.4g")

    print(f"old SNPs {len(o):,}, new SNPs {len(n):,}, matched {len(m):,}")
    print(m.allele_match.value_counts().to_string())
    print(f"flip: +1 {int((m.flip == 1).sum()):,}, -1 {int((m.flip == -1).sum()):,}, "
          f"unresolved {int(m.flip.isna().sum()):,}; allele vs genotype disagree {int((~m.agree).sum()):,}")


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:])

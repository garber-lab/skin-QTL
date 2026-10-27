"""
eQTL integration, analyses A (annotation), D (enrichment) and E (direction concordance), one cell type per call.

  python eqtl_integrate.py <cell_type> <C_dir> <align.tsv.gz> <eqtl_dir> <out_dir> [q_cutoff=0.05]

<C_dir>     output of snv_effects.py (<CT>_snv_effects.tsv.gz: one row per SNV in a peak, hg38 REF/ALT)
<align>     all.align.tsv.gz from align_alleles.py (flip = +1 same ALT, -1 swapped)
<eqtl_dir>  .../DREG/2026models: <CT>/results/result_NNN.tsv (all tested pairs) and
            results/<CT>_<COND>_eQTL_gene_BH_table.txt (eigenMT + BH per gene)

Conditions: PBS and IFNG; the eQTL of a condition is compared with the prediction of the same condition.

A. Per SNV in a peak and condition (category):
   lead          lead SNP of an eGene (q_gene < q_cutoff)
   sig           passes the eigenMT gene-level cutoff for >= 1 eGene: p_nominal <= p_thresh / Meff_gene, where
                 p_thresh = largest p_gene_eigenMT with q_gene < q_cutoff (eigenMT: p_gene = p_min * Meff)
   intermediate  tested, not sig, min nominal p <= 0.05
   tested_null   tested, min nominal p > 0.05 over all genes
   untested      not in the eQTL pairs (e.g. no gene within 500 kb)
   'sig_all' = lead + sig. Best gene = smallest nominal p; beta_aligned = beta * flip = per copy of the hg38 ALT allele
   = the allele ChromBPNet scored (eQTL dosage = count of the old VCF's ALT allele, step0.3).
D. Enrichment: predicted caQTL (identifiable, same sign in 5 folds, |pred log2| >= t) in sig_all / lead vs null SNVs;
   Fisher test, plus Cochran-Mantel-Haenszel test stratified by MAF bin x distance-to-summit bin (matching).
E. Direction concordance: sign(beta_aligned) == sign(pred log2) for SNVs with a predicted effect, binned by |pred log2|
   (same bins as the ASA validation); one-sided binomial vs 0.5; null SNVs as control.

Writes <out_dir>/<CT>_snv_annot.tsv.gz, <CT>_enrichment.tsv, <CT>_concordance_by_bin.tsv, <CT>_categories.tsv
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

CONDS = ["PBS", "IFNG"]
THRESH = [0.2, 0.3, 0.5, 1.0]
CBINS = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, np.inf]       # as asa_analyze.py
MAF_BINS = [0, 0.1, 0.2, 0.3, 0.4, 0.51]
DIST_BINS = [0, 100, 250, 500, np.inf]
NULL_P = 0.05
# result_NNN.tsv columns used (by position: snp, gene, beta_PBS, p_PBS, beta_IFNG, p_IFNG)
RES_COLS = {0: "rsid", 1: "gene", 4: "beta_PBS", 5: "p_PBS", 8: "beta_IFNG", 9: "p_IFNG"}


def cmh(tables):
    """Mantel-Haenszel OR and CMH chi2 p (continuity-corrected) over 2x2 tables [[a, b], [c, d]]."""
    num = den = dev = var = 0.0
    for (a, b), (c, d) in tables:
        n = a + b + c + d
        if n < 2:
            continue
        num += a * d / n
        den += b * c / n
        r1, r2, c1, c2 = a + b, c + d, a + c, b + d
        dev += a - r1 * c1 / n
        var += r1 * r2 * c1 * c2 / (n * n * (n - 1))
    if var == 0:
        return np.nan, np.nan
    chi2 = (abs(dev) - 0.5) ** 2 / var
    return (num / den if den > 0 else np.inf), stats.chi2.sf(chi2, 1)


def binom_p(k, n):
    return stats.binom.sf(k - 1, n, 0.5) if n else np.nan


def load_eqtl(eqtl_dir, ct, rsids):
    """all tested pairs of this cell type, restricted to the SNVs in peaks."""
    parts = []
    for f in sorted(glob.glob(os.path.join(eqtl_dir, ct, "results", "result_*.tsv"))):
        r = pd.read_csv(f, sep="\t", usecols=list(RES_COLS), dtype={0: str, 1: str})
        r.columns = [RES_COLS[i] for i in sorted(RES_COLS)]
        parts.append(r[r.rsid.isin(rsids)])
    r = pd.concat(parts, ignore_index=True)
    print(f"{ct}: {len(parts)} result files, {len(r):,} tested pairs for SNVs in peaks")
    return r


def main(ct, cdir, align, eqtl_dir, out, q_cut=0.05):
    q_cut = float(q_cut)
    os.makedirs(out, exist_ok=True)
    s = pd.read_csv(os.path.join(cdir, f"{ct}_snv_effects.tsv.gz"), sep="\t")
    al = pd.read_csv(align, sep="\t", usecols=["chrom", "pos", "rsid", "ref_new", "alt_new", "flip"])
    s = s.merge(al, on=["chrom", "pos"], how="left")
    bad = s.rsid.notna() & ((s.ref != s.ref_new) | (s.alt != s.alt_new))
    print(f"{ct}: {len(s):,} SNVs in peaks, {s.rsid.notna().sum():,} with an rsID in the alignment, "
          f"{int(bad.sum())} REF/ALT mismatches vs the phased VCF (dropped)")
    s = s[~bad].drop(columns=["ref_new", "alt_new"])
    af = s.alt_freq.clip(0, 1)
    s["maf"] = np.minimum(af, 1 - af)
    s["stratum"] = (pd.cut(s.maf, MAF_BINS, right=False, labels=False).astype(str) + "_" +
                    pd.cut(s.dist_to_summit.abs(), DIST_BINS, right=False, labels=False).astype(str))

    pairs = load_eqtl(eqtl_dir, ct, set(s.rsid.dropna()))
    annot, enr, conc, cats = [], [], [], []
    for cond in CONDS:
        bh = pd.read_csv(os.path.join(eqtl_dir, "results", f"{ct}_{cond}_eQTL_gene_BH_table.txt"), sep="\t")
        eg = bh[bh.q_gene < q_cut]
        p_thresh = eg.p_gene_eigenMT.max() if len(eg) else 0.0
        gthr = (p_thresh / eg.set_index("gene").Meff).rename("p_cut")       # per-gene nominal cutoff
        leads = set(zip(eg.lead_snp, eg.gene))
        print(f"{ct} {cond}: {len(eg):,} eGenes (q < {q_cut}), p_gene cutoff {p_thresh:.3g}")

        p = pairs[["rsid", "gene", f"beta_{cond}", f"p_{cond}"]].rename(
            columns={f"beta_{cond}": "beta", f"p_{cond}": "p"}).dropna(subset=["p"])
        p = p.join(gthr, on="gene")
        p["is_sig"] = p.p <= p.p_cut                                            # NaN p_cut (not eGene) -> False
        p["is_lead"] = [(a, b) in leads for a, b in zip(p.rsid, p.gene)]
        # best gene per SNV: a significant pair first, then the smallest p
        p = p.sort_values(["rsid", "is_sig", "p"], ascending=[True, False, True])
        best = p.drop_duplicates("rsid").set_index("rsid")
        g = p.groupby("rsid")
        best["n_genes"] = g.size()
        best["n_sig_genes"] = g.is_sig.sum()
        best["any_lead"] = g.is_lead.any()
        best["min_p"] = g.p.min()

        a = s[["cell_type", "snv_id", "rsid", "chrom", "pos", "ref", "alt", "peak_id", "dist_to_summit", "maf",
               "stratum", "identifiable", "ld_group_size", "flip", f"{cond}_beta", f"{cond}_same_sign"]].copy()
        a = a.rename(columns={f"{cond}_beta": "pred_log2", f"{cond}_same_sign": "pred_same_sign"})
        a.insert(1, "condition", cond)
        a = a.join(best[["gene", "beta", "p", "is_sig", "n_genes", "n_sig_genes", "any_lead", "min_p"]], on="rsid")
        a["beta_aligned"] = a.beta * a.flip
        a["category"] = np.select(
            [a.any_lead.fillna(False).astype(bool), a.is_sig.fillna(False).astype(bool), a.min_p <= NULL_P,
             a.min_p > NULL_P],
            ["lead", "sig", "intermediate", "tested_null"], "untested")
        a["pred_ok"] = a.identifiable & a.pred_same_sign & a.pred_log2.notna()
        annot.append(a)
        cats.append(a.groupby("category").size().rename("n_snvs").reset_index().assign(cell_type=ct, condition=cond))

        # ---------- D: enrichment of predicted caQTLs ----------
        null = a[a.category == "tested_null"]
        for test_set, sel in [("sig_all", a.category.isin(["lead", "sig"])), ("lead", a.category == "lead")]:
            x = a[sel]
            for t in THRESH:
                hit = lambda d: d.pred_ok & (d.pred_log2.abs() >= t)
                k1, n1, k0, n0 = int(hit(x).sum()), len(x), int(hit(null).sum()), len(null)
                orr, pf = stats.fisher_exact([[k1, n1 - k1], [k0, n0 - k0]], "greater") if n1 and n0 else (np.nan,) * 2
                tabs = []
                for st in set(x.stratum) & set(null.stratum):
                    xs, ns_ = x[x.stratum == st], null[null.stratum == st]
                    a1, a0 = int(hit(xs).sum()), int(hit(ns_).sum())
                    tabs.append(((a1, len(xs) - a1), (a0, len(ns_) - a0)))
                or_mh, p_mh = cmh(tabs)
                enr.append({"cell_type": ct, "condition": cond, "eqtl_set": test_set, "abs_log2_ge": t,
                            "n_eqtl": n1, "k_eqtl": k1, "frac_eqtl": k1 / n1 if n1 else np.nan,
                            "n_null": n0, "k_null": k0, "frac_null": k0 / n0 if n0 else np.nan,
                            "OR": orr, "p_fisher": pf, "OR_MH": or_mh, "p_CMH": p_mh})

        # ---------- E: direction concordance by predicted effect bin ----------
        ok = a.pred_ok & a.beta_aligned.notna() & (a.beta_aligned != 0)
        a_ok = a[ok].assign(abin=pd.cut(a[ok].pred_log2.abs(), CBINS, right=False))
        for test_set, sel in [("sig_all", a_ok.category.isin(["lead", "sig"])), ("lead", a_ok.category == "lead"),
                              ("tested_null", a_ok.category == "tested_null")]:
            x = a_ok[sel]
            same = np.sign(x.pred_log2) == np.sign(x.beta_aligned)
            for b, grp in same.groupby(x.abin, observed=False):
                k, n = int(grp.sum()), len(grp)
                conc.append({"cell_type": ct, "condition": cond, "eqtl_set": test_set, "bin": str(b),
                             "bin_lo": b.left, "n": n, "k_same": k, "concordance": k / n if n else np.nan,
                             "p_binom": binom_p(k, n)})
            for t in THRESH:                                                     # cumulative, as in Section 18
                m = x.pred_log2.abs() >= t
                k, n = int(same[m].sum()), int(m.sum())
                conc.append({"cell_type": ct, "condition": cond, "eqtl_set": test_set, "bin": f">={t}",
                             "bin_lo": t, "n": n, "k_same": k, "concordance": k / n if n else np.nan,
                             "p_binom": binom_p(k, n)})

    annot = pd.concat(annot, ignore_index=True)
    annot.to_csv(os.path.join(out, f"{ct}_snv_annot.tsv.gz"), sep="\t", index=False, float_format="%.5g")
    cats = pd.concat(cats, ignore_index=True).pivot(index=["cell_type", "category"], columns="condition",
                                                     values="n_snvs").reset_index()
    cats.to_csv(os.path.join(out, f"{ct}_categories.tsv"), sep="\t", index=False)
    enr = pd.DataFrame(enr)
    enr.to_csv(os.path.join(out, f"{ct}_enrichment.tsv"), sep="\t", index=False, float_format="%.4g")
    conc = pd.DataFrame(conc)
    conc.to_csv(os.path.join(out, f"{ct}_concordance_by_bin.tsv"), sep="\t", index=False, float_format="%.4g")

    pd.set_option("display.width", 200)
    print("\n== categories (SNVs in peaks)\n" + cats.to_string(index=False))
    print("\n== D enrichment\n" + enr.to_string(index=False))
    print("\n== E concordance (cumulative)\n" + conc[conc.bin.str.startswith(">=")].to_string(index=False))
    ex = annot[annot.rsid == "rs6844475"]
    if len(ex):
        print("\n== rs6844475\n" + ex[["condition", "snv_id", "peak_id", "pred_log2", "pred_same_sign", "gene",
                                         "beta", "flip", "beta_aligned", "p", "category"]].to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) not in (6, 7):
        sys.exit(__doc__)
    main(*sys.argv[1:])

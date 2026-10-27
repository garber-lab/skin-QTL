"""
Step 1b, part 3: compare predicted haplotype differences with allele-specific ATAC reads (runs in the container).

  python asa_analyze.py <run_dir> [--both-alleles]

--both-alleles (sensitivity analysis): drop every donor-peak that contains a het site with >= 10 reads where
  either allele has < 2 reads. Such one-allele sites are mostly genotype (imputation) errors - the donor is
  really homozygous - and show up as bands at observed log2 = +-4. Output goes to analysis_both_alleles/.

Reads <run_dir>/bam_table.tsv, counts/*.peaks.tsv, counts/*.sites.tsv and hets/<CT>/pred.tsv; writes
<run_dir>/analysis/:
  donor_peak_table.tsv.gz   one row per context x donor x peak: observed hap1/hap2 reads, predicted hap1/hap2
                            counts (same-condition and other-condition model), log2 ratios, beta-binomial p
  summary.tsv               per cell type x condition (e.g. KRT_PBS) and pooled: REF fraction, correlation,
                            direction concordance, controls
  summary_per_bam.tsv       the same per BAM (donor x cell type x condition)
  calibration.tsv           predicted vs observed effect, binned by predicted |log2 ratio|
  concordance_by_bin.tsv    direction concordance per predicted |log2| bin, fold-mean and held-out-fold prediction
  top_sites.tsv             strongest observed imbalances with predictions
  rs6844475.tsv             the pilot locus (chr4:25020562) in its heterozygous training donors
  fig_scatter.png, fig_calibration.png, fig_concordance.png, fig_concordance_by_bin.png

Held-out-chromosome check: the 6 donors' pooled reads trained the models, and each fold model was trained on all
chromosomes except its validation + test ones. The 'heldout' prediction of a donor-peak uses only the fold whose
TEST chromosomes include the peak's chromosome (never seen in training); 'train4' = mean of the other 4 folds on
the same donor-peaks. Test chromosomes per fold: TEST_CHROMS below (ChromBPNet standard human splits), or
<run_dir>/fold_test_chroms.tsv if present (2 columns: fold number, comma-separated chromosomes).
Needs pred.tsv with per-fold columns (asa_hets.py from 2026-10-27 on).
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import betaln, gammaln

MIN_READS_COR = 10        # donor-peaks used for the correlation
MIN_READS = 20            # ... for direction concordance / calibration
MIN_PRED = 0.2            # |predicted log2(hap1/hap2)| counted as a real predicted effect (~15%)
FLAT = 0.05               # |predicted log2| below this = flat peak (control)
FILTER_MIN_DEPTH = 10     # --both-alleles: judge a het site only if it has >= this many reads ...
FILTER_MIN_ALLELE = 2     # ... and call it one-allele if REF or ALT has fewer reads than this
PILOT = ("chr4", 25020562)
CBINS = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, np.inf]     # predicted |log2| bins for concordance_by_bin.tsv
TEST_CHROMS = {0: "chr1 chr3 chr6", 1: "chr2 chr8 chr9 chr16", 2: "chr4 chr11 chr12 chr15 chrY",
               3: "chr5 chr10 chr14 chr18 chr20 chr22", 4: "chr7 chr13 chr17 chr19 chr21 chrX"}
OTHER = {"PBS": "IFNG", "IFNG": "PBS"}


def wcorr(x, y, w):
    mx, my = np.average(x, weights=w), np.average(y, weights=w)
    cov = np.average((x - mx) * (y - my), weights=w)
    return cov / np.sqrt(np.average((x - mx) ** 2, weights=w) * np.average((y - my) ** 2, weights=w))


def concordance(pred, obs):
    """fraction of sites where predicted and observed log2 ratios have the same sign; one-sided binomial p."""
    keep = obs != 0
    k, n = int((np.sign(pred[keep]) == np.sign(obs[keep])).sum()), int(keep.sum())
    p = stats.binom.sf(k - 1, n, 0.5) if n else np.nan
    return k, n, (k / n if n else np.nan), p


def wilson(k, n, z=1.96):
    """95% Wilson interval of a proportion k / n."""
    if n == 0:
        return np.nan, np.nan
    p = k / n
    m = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    return c - m, c + m


def test_chroms(run):
    f = os.path.join(run, "fold_test_chroms.tsv")
    if os.path.exists(f):
        sp = {int(a): b.replace(",", " ") for a, b in (l.split()[:2] for l in open(f) if l.strip())}
        print(f"fold test chromosomes from {f}")
    else:
        sp = TEST_CHROMS
        print("fold test chromosomes: ChromBPNet standard splits (TEST_CHROMS in asa_analyze.py)")
    for k in sorted(sp):
        print(f"  fold{k}: {sp[k]}")
    return {ch: k for k, v in sp.items() for ch in v.split()}


def bb_logpmf(k, n, mu, rho):
    a, b = mu * (1 - rho) / rho, (1 - mu) * (1 - rho) / rho
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1) + betaln(k + a, n - k + b) - betaln(a, b)


def bb_two_sided(k, n, rho):
    """two-sided beta-binomial p-value vs 50:50 (sum of outcomes no more likely than the observed one)."""
    out = np.empty(len(k))
    for i, (ki, ni) in enumerate(zip(k, n)):
        lp = bb_logpmf(np.arange(ni + 1), ni, 0.5, rho)
        out[i] = min(1.0, np.exp(lp[lp <= lp[ki] + 1e-9]).sum())
    return out


def fit_rho(k, n):
    """overdispersion of the read split around 50:50, by maximum likelihood on flat peaks."""
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(lambda r: -bb_logpmf(k, n, 0.5, r).sum(), bounds=(1e-4, 0.5), method="bounded")
    return res.x


def main(run, both_alleles=False):
    out = os.path.join(run, "analysis_both_alleles" if both_alleles else "analysis")
    os.makedirs(out, exist_ok=True)
    bams = pd.read_csv(os.path.join(run, "bam_table.tsv"), sep="\t")

    # ---------- observed counts + predictions ----------
    peaks = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(os.path.join(run, "counts", "*.peaks.tsv")))])
    sites = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(os.path.join(run, "counts", "*.sites.tsv")))])
    bams["ct_cond"] = bams.cell_type + "_" + bams.condition
    peaks = peaks.merge(bams[["context", "cell_type", "condition", "ct_cond"]], on="context")
    sites = sites.merge(bams[["context", "cell_type", "condition", "ct_cond"]], on="context")
    pred = pd.concat([pd.read_csv(os.path.join(run, "hets", ct, "pred.tsv"), sep="\t").assign(cell_type=ct)
                      for ct in sorted(bams.cell_type.unique())])
    d = peaks.merge(pred, on=["cell_type", "donor", "peak_id"], how="inner")
    if len(d) != len(peaks):
        print(f"WARNING {len(peaks) - len(d)} counted peaks without predictions (dropped)")
    pc = sites.assign(peak_id=sites.peak_ids.str.split(",")).explode("peak_id")
    d = d.merge(pc[["cell_type", "peak_id", "chrom"]].drop_duplicates(["cell_type", "peak_id"]),
                on=["cell_type", "peak_id"], how="left")

    cond = d.condition.values
    pick = lambda c, h: np.where(cond == "PBS", d[f"PBS_{h}"], d[f"IFNG_{h}"]) if c == "same" else \
        np.where(cond == "PBS", d[f"IFNG_{h}"], d[f"PBS_{h}"])
    d["pred_hap1"], d["pred_hap2"] = pick("same", "hap1"), pick("same", "hap2")
    d["pred_other_hap1"], d["pred_other_hap2"] = pick("other", "hap1"), pick("other", "hap2")
    d["pred_log2"] = np.log2(d.pred_hap1 / d.pred_hap2)
    d["pred_other_log2"] = np.log2(d.pred_other_hap1 / d.pred_other_hap2)
    d["pred_share"] = d.pred_hap1 / (d.pred_hap1 + d.pred_hap2)
    d["reads"] = d.hap1_reads + d.hap2_reads
    d["obs_log2"] = np.log2((d.hap1_reads + 0.5) / (d.hap2_reads + 0.5))
    d["obs_share"] = d.hap1_reads / d.reads.where(d.reads > 0)

    # held-out-chromosome prediction (fold whose test set contains the peak's chromosome) and the other 4 folds
    folds = sorted({int(m.group(1)) for c in d.columns for m in [re.match(r"PBS_fold(\d+)_hap1$", c)] if m})
    has_folds = len(folds) > 1
    d["heldout_fold"] = d.chrom.map(test_chroms(run))
    h = {k: np.full(len(d), np.nan) for k in ("h1", "h2", "t1", "t2")}
    if has_folds:
        for f in folds:
            others = [g for g in folds if g != f]
            for c in ("PBS", "IFNG"):
                m = ((d.heldout_fold == f) & (d.condition == c)).values
                h["h1"][m], h["h2"][m] = d.loc[m, f"{c}_fold{f}_hap1"], d.loc[m, f"{c}_fold{f}_hap2"]
                h["t1"][m] = d.loc[m, [f"{c}_fold{g}_hap1" for g in others]].mean(axis=1)
                h["t2"][m] = d.loc[m, [f"{c}_fold{g}_hap2" for g in others]].mean(axis=1)
        print(f"held-out fold assigned for {int(d.heldout_fold.notna().sum())} of {len(d)} donor-peaks")
    else:
        print("NOTE pred.tsv has no per-fold columns (made by an older asa_hets.py) -> held-out check skipped; "
              "run `bash asa_step1b.sh folds`")
    d["pred_heldout_log2"] = np.log2(h["h1"] / h["h2"])
    d["pred_train4_log2"] = np.log2(h["t1"] / h["t2"])

    if both_alleles:
        # one-allele sites -> drop the whole donor-peak (its haplotype read counts include that site)
        depth = sites.ref_reads + sites.alt_reads
        bad = (depth >= FILTER_MIN_DEPTH) & ((sites.ref_reads < FILTER_MIN_ALLELE) | (sites.alt_reads < FILTER_MIN_ALLELE))
        b = sites[bad].assign(peak_id=sites.loc[bad, "peak_ids"].str.split(",")).explode("peak_id")
        drop = set(zip(b.context, b.peak_id))
        is_bad = np.array([k in drop for k in zip(d.context, d.peak_id)], dtype=bool)
        print(f"--both-alleles: {int(bad.sum())} of {int((depth >= FILTER_MIN_DEPTH).sum())} het sites with >= "
              f"{FILTER_MIN_DEPTH} reads have < {FILTER_MIN_ALLELE} reads of one allele; dropping {int(is_bad.sum())} "
              f"of {len(d)} donor-peaks ({int((is_bad & (d.reads >= MIN_READS)).sum())} of "
              f"{int((d.reads >= MIN_READS).sum())} with >= {MIN_READS} reads)")
        d = d[~is_bad].copy()
        sites = sites[~bad]

    # beta-binomial p per donor-peak; overdispersion fitted on flat peaks with enough reads
    flat = d[(d.pred_log2.abs() < FLAT) & (d.reads >= MIN_READS)]
    rho = fit_rho(flat.hap1_reads.values, flat.reads.values) if len(flat) > 50 else 0.01
    ok = d.reads >= MIN_READS
    d["bb_p"] = np.nan
    d.loc[ok, "bb_p"] = bb_two_sided(d.loc[ok, "hap1_reads"].values.astype(int), d.loc[ok, "reads"].values.astype(int), rho)
    d.drop(columns=[c for c in d.columns if re.search(r"_fold\d+_hap[12]$", c)]).to_csv(
        os.path.join(out, "donor_peak_table.tsv.gz"), sep="\t", index=False, float_format="%.4g")

    # ---------- summary per cell type x condition (and per BAM) + pooled ----------
    def summarize(key, seed=1):
        rng = np.random.default_rng(seed)
        rows = []
        for name, g in [(c, g) for c, g in d.groupby(key)] + [("ALL", d)]:
            s = sites[sites[key] == name] if name != "ALL" else sites
            rf = s.ref_reads.sum() / max(s.ref_reads.sum() + s.alt_reads.sum(), 1)
            gc = g[g.reads >= MIN_READS_COR]
            r_w = wcorr(gc.pred_log2.values, gc.obs_log2.values, gc.reads.values) if len(gc) > 2 else np.nan
            rho_s = stats.spearmanr(gc.pred_log2, gc.obs_log2)[0] if len(gc) > 2 else np.nan
            gs = g[(g.reads >= MIN_READS) & (g.pred_log2.abs() >= MIN_PRED)]
            k, n, frac, p = concordance(gs.pred_log2.values, gs.obs_log2.values)
            go = g[(g.reads >= MIN_READS) & (g.pred_other_log2.abs() >= MIN_PRED)]
            ko, no, fo, po = concordance(go.pred_other_log2.values, go.obs_log2.values)
            flip = rng.choice([-1, 1], size=len(gs))           # shuffled hap1/hap2 labels
            ks, ns, fs, ps = concordance(gs.pred_log2.values, gs.obs_log2.values * flip)
            gf = g[(g.reads >= MIN_READS) & (g.pred_log2.abs() < FLAT)]
            gh = g[(g.reads >= MIN_READS) & (g.pred_heldout_log2.abs() >= MIN_PRED)]     # held-out fold
            kh, nh, fh, ph = concordance(gh.pred_heldout_log2.values, gh.obs_log2.values)
            kt, nt, ft, pt = concordance(gh.pred_train4_log2.values, gh.obs_log2.values)  # same rows, 4 other folds
            gch = gc[gc.pred_heldout_log2.notna()]
            r_h = wcorr(gch.pred_heldout_log2.values, gch.obs_log2.values, gch.reads.values) if len(gch) > 2 else np.nan
            rows.append({key: name, "donors": g.donor.nunique(), "donor_peaks": len(g),
                         "donor_peaks_ge10reads": len(gc), "donor_peaks_ge20reads": int((g.reads >= MIN_READS).sum()),
                         "ref_fraction_sites": round(rf, 4), "weighted_r": round(r_w, 3), "spearman": round(rho_s, 3),
                         "n_concord_test": n, "concordance": round(frac, 3), "concordance_p": p,
                         "cross_condition_n": no, "cross_condition_concordance": round(fo, 3),
                         "shuffled_n": ns, "shuffled_concordance": round(fs, 3),
                         "flat_n": len(gf), "flat_pooled_hap1_share": round(gf.hap1_reads.sum() / max(gf.reads.sum(), 1), 3),
                         "overdispersion_rho": round(rho, 4),
                         "heldout_n": nh, "heldout_concordance": round(fh, 3),
                         "train4_concordance_same_rows": round(ft, 3), "heldout_weighted_r": round(r_h, 3)})
        return pd.DataFrame(rows)

    summ = summarize("ct_cond")
    summ.to_csv(os.path.join(out, "summary.tsv"), sep="\t", index=False)
    summarize("context").to_csv(os.path.join(out, "summary_per_bam.tsv"), sep="\t", index=False)

    # ---------- direction concordance per predicted-effect bin ----------
    def by_bin(col, label):
        e = d[(d.reads >= MIN_READS) & d[col].notna()]
        e = e.assign(bin=pd.cut(e[col].abs(), CBINS, right=False))
        rows = []
        for (ctx, b), g in list(e.groupby(["ct_cond", "bin"], observed=True)) + \
                [(("ALL", b), g) for b, g in e.groupby("bin", observed=True)]:
            k, n, frac, p = concordance(g[col].values, g.obs_log2.values)
            lo, hi = wilson(k, n)
            rows.append(dict(context=ctx, prediction=label, pred_bin=str(b), n=n, concordance=frac,
                             ci_low=lo, ci_high=hi, p=p))
        return rows
    cb = pd.DataFrame(by_bin("pred_log2", "fold_mean") + (by_bin("pred_heldout_log2", "heldout_fold") if has_folds else []))
    cb.to_csv(os.path.join(out, "concordance_by_bin.tsv"), sep="\t", index=False, float_format="%.4g")

    # ---------- calibration: orient each donor-peak so the predicted-higher haplotype is 'A' ----------
    c = d[d.reads >= MIN_READS].copy()
    sgn = np.sign(c.pred_log2).replace(0, 1)
    c["pred_abs"] = c.pred_log2.abs()
    c["reads_A"] = np.where(sgn > 0, c.hap1_reads, c.hap2_reads)
    c["reads_B"] = c.reads - c.reads_A
    bins = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, np.inf]
    c["bin"] = pd.cut(c.pred_abs, bins, right=False)
    cal = []
    for (ctx, b), g in list(c.groupby(["ct_cond", "bin"], observed=True)) + \
            [(("ALL", b), g) for b, g in c.groupby("bin", observed=True)]:
        cal.append(dict(context=ctx, pred_bin=str(b), n=len(g), mean_pred_log2=g.pred_abs.mean(),
                        pooled_obs_log2=np.log2((g.reads_A.sum() + 0.5) / (g.reads_B.sum() + 0.5)),
                        obs_share_A=g.reads_A.sum() / g.reads.sum()))
    cal = pd.DataFrame(cal)
    # one-number calibration: weighted slope through the origin of observed on predicted log2 ratio
    w = c.reads.values
    slope = float((w * c.pred_log2 * c.obs_log2).sum() / (w * c.pred_log2 ** 2).sum())
    ch = c[c.pred_heldout_log2.notna()]
    slope_h = float((ch.reads * ch.pred_heldout_log2 * ch.obs_log2).sum() / (ch.reads * ch.pred_heldout_log2 ** 2).sum()) \
        if len(ch) else np.nan
    cal.to_csv(os.path.join(out, "calibration.tsv"), sep="\t", index=False, float_format="%.4g")

    # ---------- top sites and the pilot locus ----------
    cols = ["context", "donor", "peak_id", "n_het_core", "hap1_reads", "hap2_reads", "obs_log2",
            "pred_log2", "pred_other_log2", "bb_p"]
    d[d.bb_p.notna()].sort_values("bb_p")[cols].head(200).to_csv(os.path.join(out, "top_sites.tsv"), sep="\t",
                                                                   index=False, float_format="%.4g")
    pil = sites[(sites.chrom == PILOT[0]) & (sites.pos == PILOT[1])].copy()
    pil["peak_id"] = pil.peak_ids.str.split(",")
    pil = pil.explode("peak_id")[["context", "donor", "peak_id", "ref", "alt", "hap1_alt", "ref_reads", "alt_reads"]]
    pil_peaks = pil.merge(d[cols], on=["context", "donor", "peak_id"], how="left")
    pil_peaks.to_csv(os.path.join(out, "rs6844475.tsv"), sep="\t", index=False, float_format="%.4g")

    # ---------- figures ----------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ctxs = sorted(d.ct_cond.unique())
    ncol = 3
    nrow = int(np.ceil(len(ctxs) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.8 * nrow), squeeze=False)
    for ax, ctx in zip(axs.ravel(), ctxs):
        g = d[(d.ct_cond == ctx) & (d.reads >= MIN_READS)]
        ax.hexbin(g.pred_log2.clip(-2, 2), g.obs_log2.clip(-4, 4), gridsize=40, bins="log", mincnt=1, cmap="Blues")
        ax.axhline(0, lw=0.5, c="grey"); ax.axvline(0, lw=0.5, c="grey")
        s = summ.set_index("ct_cond").loc[ctx]
        ax.set_title(f"{ctx} ({int(s.donors)} donors)\nr_w = {s.weighted_r}, concordance {s.concordance}", fontsize=9)
        ax.set_xlabel("predicted log2(hap1/hap2)"); ax.set_ylabel(f"observed log2 (≥{MIN_READS} reads)")
    for ax in axs.ravel()[len(ctxs):]:
        ax.axis("off")
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig_scatter.png"), dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 4.5))
    for ctx, g in cal.groupby("context"):
        ax.plot(g.mean_pred_log2, g.pooled_obs_log2, "o-" if ctx == "ALL" else ".-",
                lw=2.5 if ctx == "ALL" else 1, c="black" if ctx == "ALL" else None, label=ctx, alpha=1 if ctx == "ALL" else 0.6)
    lim = max(1.0, cal.mean_pred_log2.max())
    ax.plot([0, lim], [0, lim], "--", c="grey", lw=0.8, label="y = x")
    ax.set_xlabel("predicted |log2| (higher vs lower haplotype)"); ax.set_ylabel("observed log2, pooled reads")
    ax.set_title(f"Calibration (slope {slope:.2f})", fontsize=10); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig_calibration.png"), dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(1.2 * len(summ) + 2, 4))
    x = np.arange(len(summ))
    for i, (col, lab) in enumerate([("concordance", "matched model"), ("cross_condition_concordance", "other condition"),
                                    ("shuffled_concordance", "shuffled haplotypes")]):
        ax.bar(x + (i - 1) * 0.27, summ[col], 0.27, label=lab)
    ax.axhline(0.5, ls="--", c="grey", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels(summ.ct_cond, rotation=30, ha="right", fontsize=8)
    ax.set_ylim(0.3, 1); ax.set_ylabel(f"direction concordance\n(|pred log2| ≥ {MIN_PRED}, ≥ {MIN_READS} reads)")
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(os.path.join(out, "fig_concordance.png"), dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4.2))
    labs = [str(b) for b in pd.cut([], CBINS, right=False).categories]
    for (ctx, lab), g in cb.groupby(["context", "prediction"]):
        xi = [labs.index(b) for b in g.pred_bin]
        if ctx == "ALL":
            ax.errorbar(np.array(xi) + (0.08 if lab == "heldout_fold" else 0), g.concordance,
                        yerr=[g.concordance - g.ci_low, g.ci_high - g.concordance], fmt="o-" if lab == "fold_mean" else "s--",
                        c="black" if lab == "fold_mean" else "tab:orange", lw=2, capsize=3,
                        label="all contexts, mean of 5 folds" if lab == "fold_mean" else "all contexts, held-out fold only")
        elif lab == "fold_mean":
            ax.plot(xi, g.concordance, ".-", lw=0.8, alpha=0.5, label=ctx)
    ax.axhline(0.5, ls="--", c="grey", lw=0.8)
    ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, rotation=30, ha="right", fontsize=8)
    ax.set_xlabel("predicted |log2(hap1/hap2)|"); ax.set_ylabel(f"direction concordance (≥ {MIN_READS} reads)")
    ax.legend(fontsize=6.5); fig.tight_layout(); fig.savefig(os.path.join(out, "fig_concordance_by_bin.png"), dpi=150)
    plt.close(fig)

    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(summ.to_string(index=False))
    print(f"\nCalibration slope (observed log2 per predicted log2, all contexts): {slope:.3f}"
          f"; held-out fold only: {slope_h:.3f}")
    print("\nConcordance per predicted |log2| bin, all contexts:")
    print(cb[cb.context == "ALL"].to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print(f"Pilot locus rows: {len(pil_peaks)}; outputs in {out}")


if __name__ == "__main__":
    main(sys.argv[1], "--both-alleles" in sys.argv[2:])

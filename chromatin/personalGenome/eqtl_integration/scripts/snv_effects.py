"""
eQTL integration, analysis C: per-SNV predicted accessibility effect (predicted caQTLs) from the haplotype predictions.

  python snv_effects.py <merged_dir> <cell_type> <out_dir>

<merged_dir> = .../multipeak/<CT>_donors_phasedMAF0.05/merged (haplotypes.tsv.gz, peak_info.tsv).

Method (per peak): the peak's unique donor haplotypes differ by sets of ALT alleles (hg38 REF/ALT, alt_alleles column
'pos REF>ALT', 1-based VCF POS). Regress log2(predicted count) of each haplotype on 0/1 ALT indicators of the peak's
SNVs (+ intercept), weighted by the number of donor haplotypes carrying it (weight >= 1, so the REF haplotype counts
even if no donor carries it). Coefficient = additive log2 effect of ALT vs REF on that background.

  log2(count_h) = a + sum_j beta_j * x_hj        x_hj = 1 if haplotype h carries the ALT of SNV j

- beta      : fit on log2 of the fold-mean count (*_mean column) = same scale as the ASA validation (Section 18)
- beta_foldK: fit on each fold separately; same_sign = all 5 fold betas have the sign of beta
- SNVs whose ALT always occurs together on the observed haplotypes (identical indicator columns) cannot be separated:
  they share one coefficient (ld_group_size > 1, ld_group lists them). Columns that are otherwise not identifiable
  (e.g. x_A = x_B + x_C) get beta = NaN, identifiable = False.
- pair_beta : model-free check: mean log2 ratio over haplotype pairs that differ ONLY at this SNV (group); n_pairs.
- max_resid : largest |residual| of the additive fit in the peak (PBS_mean); large = non-additive (epistatic) peak.

Thresholds (Section 18, ASA concordance by bin, fold mean): |log2| >= 0.2 / 0.3 / 0.5 / 1.0 -> 0.73 / 0.81 / 0.88 /
0.97 of read-based directions agree. Main set: |beta| >= 0.3 and same_sign.

Writes <out_dir>/:
  <CT>_snv_effects_all.tsv.gz  one row per peak x SNV
  <CT>_snv_effects.tsv.gz      one row per SNV (the peak whose summit is closest; n_peaks = peaks containing it)
  <CT>_snv_summary.tsv         number of SNVs above each threshold, per condition
"""
import os
import re
import sys

import numpy as np
import pandas as pd

CONDS = ["PBS", "IFNG"]
FOLDS = range(5)
THRESH = [0.2, 0.3, 0.5, 1.0]
PRECISION = {0.2: 0.73, 0.3: 0.81, 0.5: 0.88, 1.0: 0.97}   # ASA concordance at |pred log2| >= t (Section 18)
MAIN = 0.3
SNV_RE = re.compile(r"(\d+)([ACGT]+)>([ACGT]+)")


def parse_alts(s):
    """'778639A>G;779047G>A' -> [(778639,'A','G'), ...]; 'REF' -> []. Any separator between entries is accepted."""
    return [(int(p), r, a) for p, r, a in SNV_RE.findall(s)] if isinstance(s, str) else []


def wls(Z, Y, w):
    """weighted least squares, columns of Y fitted separately; returns coefficients (ncol x nY) and residuals."""
    sw = np.sqrt(w)[:, None]
    coef, *_ = np.linalg.lstsq(Z * sw, Y * sw, rcond=None)
    return coef, Y - Z @ coef


def identifiable_cols(Z):
    """True for columns whose coefficient is uniquely determined (dropping the column lowers the rank)."""
    r = np.linalg.matrix_rank(Z)
    if r == Z.shape[1]:
        return np.ones(Z.shape[1], bool)
    return np.array([np.linalg.matrix_rank(np.delete(Z, j, axis=1)) < r for j in range(Z.shape[1])])


def fit_peak(peak_id, hap_snvs, n, Y):
    """
    hap_snvs: list (one per haplotype) of SNV tuples; n: donor haplotypes per haplotype; Y: haplotypes x outcomes
    (log2 counts, columns = OUT_COLS). Returns one dict per SNV.
    """
    snvs = sorted({s for hs in hap_snvs for s in hs})
    if not snvs:
        return []
    col = {s: j for j, s in enumerate(snvs)}
    X = np.zeros((len(hap_snvs), len(snvs)), np.int8)
    for h, hs in enumerate(hap_snvs):
        for s in hs:
            X[h, col[s]] = 1

    # SNVs with identical indicator columns = inseparable on these haplotypes -> one shared coefficient
    uniq, grp = np.unique(X.T, axis=0, return_inverse=True)
    grp = grp.ravel()
    U = uniq.T.astype(float)                                   # haplotypes x groups
    Z = np.column_stack([np.ones(len(U)), U])
    w = np.maximum(n, 1).astype(float)
    ident = identifiable_cols(Z)[1:]
    coef, resid = wls(Z, Y, w)
    beta = coef[1:]                                            # groups x outcomes
    beta[~ident] = np.nan
    max_resid = float(np.abs(resid[:, OUT_COLS.index("PBS_mean")]).max())

    # model-free check: haplotype pairs that differ only in group g
    pair = np.full((U.shape[1], len(OUT_COLS)), np.nan)
    npairs = np.zeros(U.shape[1], int)
    for g in range(U.shape[1]):
        rest = np.delete(U, g, axis=1)
        keys = {}
        for h in range(len(U)):
            keys.setdefault(rest[h].tobytes(), []).append(h)
        diffs = []
        for hs in keys.values():
            on = [h for h in hs if U[h, g] == 1]
            off = [h for h in hs if U[h, g] == 0]
            diffs += [Y[a] - Y[b] for a in on for b in off]
        if diffs:
            pair[g] = np.mean(diffs, axis=0)
            npairs[g] = len(diffs)

    total = n.sum()
    rows = []
    for j, s in enumerate(snvs):
        g = grp[j]
        members = [snvs[k] for k in np.where(grp == g)[0]]
        r = {"peak_id": peak_id, "pos": s[0], "ref": s[1], "alt": s[2],
             "n_haps_peak": len(U), "n_snvs_peak": len(snvs),
             "alt_hap_count": int(n[X[:, j] == 1].sum()), "alt_freq": n[X[:, j] == 1].sum() / total,
             "ld_group_size": len(members),
             "ld_group": ",".join(f"{p}{a}>{b}" for p, a, b in members) if len(members) > 1 else "",
             "identifiable": bool(ident[g]), "n_pairs": npairs[g], "max_resid": max_resid}
        for k, c in enumerate(OUT_COLS):
            r[c] = beta[g, k]
        for cond in CONDS:
            r[f"{cond}_pair_beta"] = pair[g, OUT_COLS.index(f"{cond}_mean")]
        rows.append(r)
    return rows


OUT_COLS = [f"{c}_fold{k}" for c in CONDS for k in FOLDS] + [f"{c}_mean" for c in CONDS]


def main(merged, ct, out):
    os.makedirs(out, exist_ok=True)
    hap = pd.read_csv(os.path.join(merged, "haplotypes.tsv.gz"), sep="\t")
    info = pd.read_csv(os.path.join(merged, "peak_info.tsv"), sep="\t", usecols=["peak_id", "chrom", "center"])
    print(f"{ct}: {len(hap):,} haplotype rows, {hap.peak_id.nunique():,} peaks")

    hap = hap.sort_values(["peak_id", "hap_id"], kind="stable").reset_index(drop=True)
    snv_lists = [parse_alts(s) for s in hap.alt_alleles]
    Yall = np.log2(hap[OUT_COLS].to_numpy(float))
    nall = hap.n_donor_haplotypes.to_numpy()
    pid = hap.peak_id.to_numpy()
    starts = np.r_[0, np.where(pid[1:] != pid[:-1])[0] + 1, len(pid)]

    rows = []
    for i, (a, b) in enumerate(zip(starts[:-1], starts[1:])):
        rows += fit_peak(pid[a], snv_lists[a:b], nall[a:b], Yall[a:b])
        if i % 20000 == 0:
            print(f"  {i:,} peaks", flush=True)
    d = pd.DataFrame(rows).rename(columns={f"{c}_mean": f"{c}_beta" for c in CONDS})
    d = d.rename(columns={f"{c}_fold{k}": f"{c}_beta_fold{k}" for c in CONDS for k in FOLDS})

    d = d.merge(info, on="peak_id", how="left")
    d.insert(0, "cell_type", ct)
    d["snv_id"] = d.chrom + ":" + d.pos.astype(str) + ":" + d.ref + ":" + d.alt
    d["dist_to_summit"] = (d.pos - 1) - d.center                      # center = 0-based summit
    d["in_central_1kb"] = d.dist_to_summit.abs() < 500
    for c in CONDS:
        folds = d[[f"{c}_beta_fold{k}" for k in FOLDS]].to_numpy()
        d[f"{c}_same_sign"] = (np.sign(folds) == np.sign(d[f"{c}_beta"].to_numpy())[:, None]).all(axis=1)
        ok = d[f"{c}_same_sign"] & d[f"{c}_beta"].notna()
        d[f"{c}_precision"] = np.nan
        for t in THRESH:                                                  # ascending -> highest bin wins
            d.loc[ok & (d[f"{c}_beta"].abs() >= t), f"{c}_precision"] = PRECISION[t]
        d[f"{c}_caqtl"] = ok & (d[f"{c}_beta"].abs() >= MAIN)
    d["IFNG_minus_PBS_beta"] = d.IFNG_beta - d.PBS_beta
    d["caqtl_any"] = d.PBS_caqtl | d.IFNG_caqtl

    lead = ["cell_type", "snv_id", "chrom", "pos", "ref", "alt", "peak_id", "dist_to_summit", "in_central_1kb",
            "alt_freq", "PBS_beta", "IFNG_beta", "IFNG_minus_PBS_beta", "PBS_same_sign", "IFNG_same_sign",
            "PBS_precision", "IFNG_precision", "PBS_caqtl", "IFNG_caqtl", "caqtl_any", "identifiable",
            "ld_group_size", "ld_group", "PBS_pair_beta", "IFNG_pair_beta", "n_pairs", "max_resid"]
    d = d[lead + [c for c in d.columns if c not in lead and c != "center"]]
    d.to_csv(os.path.join(out, f"{ct}_snv_effects_all.tsv.gz"), sep="\t", index=False, float_format="%.5g")

    # one row per SNV: prefer identifiable fits, then the peak whose summit is closest
    d["_absd"] = d.dist_to_summit.abs()
    d["_noid"] = ~d.identifiable
    n_peaks = d.groupby("snv_id").size()
    s = d.sort_values(["snv_id", "_noid", "_absd"]).drop_duplicates("snv_id").drop(columns=["_absd", "_noid"])
    s.insert(s.columns.get_loc("peak_id") + 1, "n_peaks", s.snv_id.map(n_peaks).values)
    s.to_csv(os.path.join(out, f"{ct}_snv_effects.tsv.gz"), sep="\t", index=False, float_format="%.5g")

    summ = []
    for c in CONDS:
        b, ok = s[f"{c}_beta"], s[f"{c}_same_sign"] & s[f"{c}_beta"].notna()
        for t in THRESH:
            k = int((ok & (b.abs() >= t)).sum())
            summ.append({"cell_type": ct, "condition": c, "abs_log2_ge": t, "fold_change": round(2 ** t, 2),
                         "expected_precision": PRECISION[t], "n_snvs": k, "n_tested": len(s),
                         "frac": k / len(s), "n_identifiable": int(s.identifiable.sum()),
                         "n_in_ld_group": int((s.ld_group_size > 1).sum())})
    summ = pd.DataFrame(summ)
    summ.to_csv(os.path.join(out, f"{ct}_snv_summary.tsv"), sep="\t", index=False, float_format="%.4g")
    print(summ.to_string(index=False))
    r = s[["PBS_beta", "PBS_pair_beta"]].dropna()
    print(f"regression vs pair estimate (PBS): r = {np.corrcoef(r.PBS_beta, r.PBS_pair_beta)[0, 1]:.3f} "
          f"on {len(r):,} SNVs with a clean pair")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])

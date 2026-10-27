"""
eQTL integration, analysis B: landscape of predicted genetic variability in accessibility across the 101 own donors.

  python landscape.py <multipeak_dir> <out_dir> [CT ...]

<multipeak_dir> = .../personalGenome/multipeak (reads <CT>_donors_phasedMAF0.05/merged/diploid_<COND>_mean.tsv.gz and
haplotypes.tsv.gz). Cell types default to KRT FRB MEL.

Per peak and condition (PBS, IFNG):
  donor_sd     SD across donors of log2(predicted diploid count)   (diploid = mean of the 2 haplotypes, fold mean)
  donor_range  max - min across donors of log2(predicted count)    = largest difference between two donors
  hap_max_eff  max |log2(hap / REF hap)| over haplotypes carried by >= 1 donor = largest single-haplotype effect
Fraction of peaks with donor_range / hap_max_eff >= 0.2 / 0.3 / 0.5 / 1.0 (thresholds from Section 18).

Writes <out_dir>/landscape_per_peak.tsv.gz, landscape_summary.tsv, fig_landscape.png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CONDS = ["PBS", "IFNG"]
THRESH = [0.2, 0.3, 0.5, 1.0]


def main(base, out, cts):
    os.makedirs(out, exist_ok=True)
    per, summ = [], []
    for ct in cts:
        merged = os.path.join(base, f"{ct}_donors_phasedMAF0.05", "merged")
        hap = pd.read_csv(os.path.join(merged, "haplotypes.tsv.gz"), sep="\t",
                          usecols=["peak_id", "hap_id", "n_donor_haplotypes", "PBS_mean", "IFNG_mean"])
        for cond in CONDS:
            m = pd.read_csv(os.path.join(merged, f"diploid_{cond}_mean.tsv.gz"), sep="\t", index_col=0)
            l2 = np.log2(m)
            d = pd.DataFrame({"cell_type": ct, "condition": cond, "peak_id": m.index,
                              "mean_log2": l2.mean(axis=1).values, "donor_sd": l2.std(axis=1).values,
                              "donor_range": (l2.max(axis=1) - l2.min(axis=1)).values})
            # largest effect of a carried haplotype vs the REF haplotype (hap_id 0)
            ref = hap[hap.hap_id == 0].set_index("peak_id")[f"{cond}_mean"]
            h = hap[(hap.hap_id != 0) & (hap.n_donor_haplotypes > 0)]
            eff = np.log2(h[f"{cond}_mean"].values / ref.reindex(h.peak_id).values)
            d["hap_max_eff"] = d.peak_id.map(pd.Series(np.abs(eff), index=h.peak_id.values).groupby(level=0).max())
            per.append(d)
            for t in THRESH:
                summ.append({"cell_type": ct, "condition": cond, "abs_log2_ge": t, "fold_change": round(2 ** t, 2),
                             "n_peaks": len(d),
                             "frac_donor_range": (d.donor_range >= t).mean(),
                             "frac_hap_max_eff": (d.hap_max_eff >= t).mean()})
            print(f"{ct} {cond}: {len(d):,} peaks, median donor_range {d.donor_range.median():.3f}, "
                  f"median donor_sd {d.donor_sd.median():.3f}", flush=True)

    per = pd.concat(per, ignore_index=True)
    per.to_csv(os.path.join(out, "landscape_per_peak.tsv.gz"), sep="\t", index=False, float_format="%.5g")
    summ = pd.DataFrame(summ)
    summ.to_csv(os.path.join(out, "landscape_summary.tsv"), sep="\t", index=False, float_format="%.4f")
    print(summ.to_string(index=False))

    # histogram of donor_range per cell type x condition, threshold lines
    fig, axes = plt.subplots(len(cts), 2, figsize=(9, 2.6 * len(cts)), sharex=True, squeeze=False)
    bins = np.linspace(0, 1.5, 61)
    for i, ct in enumerate(cts):
        for j, cond in enumerate(CONDS):
            ax = axes[i, j]
            x = per[(per.cell_type == ct) & (per.condition == cond)].donor_range
            ax.hist(np.clip(x, 0, 1.5), bins=bins, color="#4C72B0")
            for t in THRESH[:3]:
                ax.axvline(t, color="#C44E52", lw=0.8, ls="--")
                ax.text(t, ax.get_ylim()[1] * 0.9, f" {(x >= t).mean():.1%}", fontsize=7, color="#C44E52")
            ax.set_title(f"{ct} {cond} (n = {len(x):,} peaks)", fontsize=9)
            ax.set_yscale("log")
    for ax in axes[-1]:
        ax.set_xlabel("donor range of log2 predicted accessibility (max - min; clipped at 1.5)")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_landscape.png"), dpi=150)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3:] or ["KRT", "FRB", "MEL"])

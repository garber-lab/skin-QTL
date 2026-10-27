#!/usr/bin/env python
"""
Motif co-occurrence and spacing from Fi-NeMo all-motif hit tables (step7_motifScan_allhits.sh).
This is the hit-based half of the cooperativity question: which motif pairs occupy the same peak
more often than chance, and do they sit at preferred distances and orientations?

Usage:
  python step7_motifScan_pairs.py --hits <dir>/<model>/hits.tsv --model KRT_IFN --out <dir>/pairs
  (one call per model; or pass several --hits/--model pairs by running it repeatedly)

Outputs, under <out>/:
  <model>.motif_summary.tsv   per motif: instances, peaks containing it, mean coefficient
  <model>.pairs.tsv           per motif pair: peaks with both, expected under independence,
                              observed/expected ratio, Fisher p, and median |distance|
  <model>.spacing_<A>_<B>.tsv distance histogram for the requested pair(s), with orientation

Counting rules:
  - a peak "contains" a motif if it has >= 1 hit of it (peak_id groups the hits)
  - expected co-occurrence assumes independence given each motif's peak frequency
  - distance is between hit midpoints, signed by the order along the genome, and orientation is
    same-strand or opposite-strand (from Fi-NeMo's strand column)
"""
import argparse, itertools, os, sys
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

p = argparse.ArgumentParser()
p.add_argument("--hits", required=True)
p.add_argument("--model", required=True)
p.add_argument("--out", required=True)
p.add_argument("--min_peaks", type=int, default=50, help="skip motifs present in fewer peaks")
p.add_argument("--min_similarity", type=float, default=0.0, help="drop weak hits below this cosine similarity")
p.add_argument("--top_pairs", type=int, default=20, help="write spacing tables for this many top pairs")
p.add_argument("--spacing_max", type=int, default=300, help="max |distance| (bp) kept in spacing tables")
p.add_argument("--pairs", default=None, help="comma-separated A:B pairs to force into the spacing output")
p.add_argument("--include_neg", action="store_true", help="also use negative patterns (default: positive only)")
a = p.parse_args()

os.makedirs(a.out, exist_ok=True)
h = pd.read_csv(a.hits, sep="\t")
need = {"peak_id", "motif_name", "start", "end", "strand", "hit_coefficient", "hit_similarity"}
missing = need - set(h.columns)
if missing:
    sys.exit(f"ERROR: {a.hits} is missing columns: {sorted(missing)}")

if a.min_similarity > 0:
    h = h[h["hit_similarity"] >= a.min_similarity]
if not a.include_neg:
    h = h[~h["motif_name"].astype(str).str.endswith("_neg")]
if h.empty:
    sys.exit("ERROR: no hits left after filtering")
h["mid"] = (h["start"] + h["end"]) / 2.0
n_peaks = h["peak_id"].nunique()
print(f"{a.model}: {len(h)} hits, {h['motif_name'].nunique()} motifs, {n_peaks} peaks with >=1 hit")

# ---------------- per-motif summary ----------------
summ = (h.groupby("motif_name")
         .agg(instances=("motif_name", "size"), peaks=("peak_id", "nunique"),
              mean_coefficient=("hit_coefficient", "mean"), mean_similarity=("hit_similarity", "mean"))
         .reset_index().sort_values("instances", ascending=False))
summ.insert(0, "model", a.model)
summ["peak_frequency"] = summ["peaks"] / n_peaks
summ.to_csv(os.path.join(a.out, f"{a.model}.motif_summary.tsv"), sep="\t", index=False)

keep = summ[summ["peaks"] >= a.min_peaks]["motif_name"].tolist()
print(f"{len(keep)} motifs in >= {a.min_peaks} peaks")
if len(keep) < 2:
    sys.exit("ERROR: fewer than 2 motifs pass --min_peaks")

# ---------------- pairwise co-occurrence ----------------
peaks_of = {mth: set(g["peak_id"]) for mth, g in h[h["motif_name"].isin(keep)].groupby("motif_name")}
rows = []
for A, B in itertools.combinations(sorted(keep), 2):
    pa, pb = peaks_of[A], peaks_of[B]
    both = len(pa & pb)
    only_a, only_b = len(pa) - both, len(pb) - both
    neither = n_peaks - both - only_a - only_b
    expected = len(pa) * len(pb) / n_peaks
    odds, pval = fisher_exact([[both, only_a], [only_b, neither]], alternative="greater")
    rows.append(dict(model=a.model, motif_A=A, motif_B=B, peaks_A=len(pa), peaks_B=len(pb),
                     peaks_both=both, expected_both=expected,
                     obs_over_exp=(both / expected if expected > 0 else np.nan),
                     odds_ratio=odds, fisher_p=pval))
pairs = pd.DataFrame(rows)
# Benjamini-Hochberg across all tested pairs
pairs = pairs.sort_values("fisher_p").reset_index(drop=True)
pairs["fdr"] = (pairs["fisher_p"] * len(pairs) / (pairs.index + 1)).clip(upper=1.0)[::-1].cummin()[::-1]

# ---------------- spacing for co-occurring pairs ----------------
def spacing(A, B):
    """signed distance between hit midpoints of A and B within the same peak, + orientation"""
    ha = h[h["motif_name"] == A][["peak_id", "mid", "strand"]]
    hb = h[h["motif_name"] == B][["peak_id", "mid", "strand"]]
    mg = ha.merge(hb, on="peak_id", suffixes=("_A", "_B"))
    if mg.empty:
        return mg
    mg["distance"] = mg["mid_B"] - mg["mid_A"]
    mg = mg[mg["distance"].abs() <= a.spacing_max]
    mg["orientation"] = np.where(mg["strand_A"] == mg["strand_B"], "same", "opposite")
    return mg

med = []
for _, r in pairs.iterrows():
    mg = spacing(r.motif_A, r.motif_B)
    med.append(mg["distance"].abs().median() if len(mg) else np.nan)
pairs["median_abs_distance"] = med
pairs.to_csv(os.path.join(a.out, f"{a.model}.pairs.tsv"), sep="\t", index=False)

wanted = [tuple(x.split(":")) for x in a.pairs.split(",")] if a.pairs else []
top = [(r.motif_A, r.motif_B) for _, r in pairs.head(a.top_pairs).iterrows()]
for A, B in dict.fromkeys(wanted + top):
    mg = spacing(A, B)
    if mg.empty:
        continue
    hist = (mg.groupby(["distance", "orientation"]).size().reset_index(name="n")
              .sort_values("distance"))
    hist.insert(0, "motif_B", B); hist.insert(0, "motif_A", A); hist.insert(0, "model", a.model)
    hist.to_csv(os.path.join(a.out, f"{a.model}.spacing_{A}_{B}.tsv"), sep="\t", index=False)

print(pairs.head(15)[["motif_A", "motif_B", "peaks_both", "expected_both",
                      "obs_over_exp", "fdr", "median_abs_distance"]].to_string(index=False))
print(f"\nWrote {a.out}/{a.model}.motif_summary.tsv, .pairs.tsv and spacing tables")

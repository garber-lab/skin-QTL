#!/usr/bin/env python
"""
Summarize per-fold variant-scorer outputs for one variant into tidy TSVs.

Inputs per model x fold (from run_snp_contrib.sh):
  per_fold/<model>_fold<f>.variant_shap.counts.h5    (deepdish)
  per_fold/<model>_fold<f>.variant_shap.profile.h5   (deepdish)
  per_fold/<model>_fold<f>.variant_predictions.h5    (h5py)
  per_fold/<model>_fold<f>.variant_scores.tsv

Outputs (<outdir>/<rsid>.*):
  contrib_window.tsv  model, shap_type, allele, pos_rel, nucleotide, is_ref_base, hyp, proj
  pred_profiles.tsv   model, allele, pos_rel, pred      (fold-averaged predicted cut counts)
  bigwig_tracks.tsv   model, track, pos_rel, value, value_cpm  (raw counts and depth-normalised;
                      value_cpm = value * 1e6 / genome-wide sum of that bigWig)
  scores_by_fold.tsv  model, fold, logfc, jsd, allele1_pred_counts, allele2_pred_counts
  scores_summary.tsv  model, mean_log2fc, sd_log2fc, alt_over_ref, folds_pos, n_folds, mean_jsd
                      (+ ism_log2fc, ism_rank / ism_n substitutions, ism_pos_rank / ism_n_pos positions)
  ism.tsv             model, pos_rel, ref_base, mut_base, mean_log2fc, sd_log2fc, is_variant, rank_abs, pos_rank
"""
import argparse, os, sys, glob
import numpy as np
import pandas as pd
import h5py
import deepdish as dd
import pyBigWig
from pyfaidx import Fasta

NUC = np.array(list("ACGT"))  # chrombpnet one-hot order

p = argparse.ArgumentParser()
p.add_argument("--outdir", required=True)
p.add_argument("--rsid", required=True)
p.add_argument("--cb", required=True)
p.add_argument("--bias", default="0.5")
p.add_argument("--models", required=True)
p.add_argument("--folds", default="0,1,2,3,4")
p.add_argument("--genome", required=True)
p.add_argument("--contrib_flank", type=int, default=50)   # +/- bp for contribution scores
p.add_argument("--track_flank", type=int, default=250)    # +/- bp for profiles/bigwigs (max 500 for profiles)
p.add_argument("--pred_bw_pattern", default=None,
               help="Optional, e.g. '{cb}/ATAC_{model}_merged_alldonors/biasFactor{bias}/pred/xxx_chrombpnet_nobias.bw'")
a = p.parse_args()

models = a.models.split(",")
folds = a.folds.split(",")
pf = os.path.join(a.outdir, "per_fold")

var = pd.read_csv(os.path.join(a.outdir, f"{a.rsid}.variant.tsv"), sep="\t", header=None,
                  names=["chr", "pos", "ref", "alt", "id"]).iloc[0]
chrom, pos, ref, alt = var.chr, int(var.pos), var.ref.upper(), var.alt.upper()
snp0 = pos - 1  # 0-based genomic coordinate

# ---- sanity check: REF matches the genome ----
gbase = str(Fasta(a.genome)[chrom][snp0]).upper()
if len(ref) == 1 and gbase != ref:
    # REF/ALT swapped relative to hg38: variant-scorer substitutes both alleles into the sequence,
    # so the analysis is unaffected; the centring checks below still verify REF/ALT placement
    print(f"WARNING: genome base at {chrom}:{pos} is {gbase}, VCF REF is {ref} "
          f"(REF/ALT swapped relative to hg38)")
print(f"{a.rsid} {chrom}:{pos} {ref}>{alt} (genome base {gbase})")

contrib_rows, pred_rows, score_rows = [], [], []

for m in models:
    # ---------------- contribution scores ----------------
    for st in ["counts", "profile"]:
        hyp_folds, seq = [], None
        for f in folds:
            fn = os.path.join(pf, f"{m}_fold{f}.variant_shap.{st}.h5")
            d = dd.io.load(fn)
            alleles = np.asarray(d["alleles"])            # 0 = allele1 (REF), 1 = allele2 (ALT)
            raw = np.asarray(d["raw"]["seq"]).astype(float)       # (2, 4, L)
            hyp = np.asarray(d["shap"]["seq"]).astype(float)      # hypothetical contributions
            order = np.argsort(alleles, kind="stable")
            raw, hyp = raw[order], hyp[order]
            if seq is None:
                seq = raw
            hyp_folds.append(hyp)
        hyp = np.mean(hyp_folds, axis=0)                   # fold-averaged, (2, 4, L)
        L = seq.shape[-1]
        c = L // 2                                          # variant-scorer puts the SNP at L//2

        # sanity checks on sequence centering
        ref_seq_base = NUC[seq[0, :, c].argmax()]
        alt_seq_base = NUC[seq[1, :, c].argmax()]
        if len(ref) == 1 and len(alt) == 1:
            assert ref_seq_base == ref and alt_seq_base == alt, \
                f"{m} {st}: centre bases {ref_seq_base}/{alt_seq_base} != {ref}/{alt}"
            diff = np.where((seq[0] != seq[1]).any(axis=0))[0]
            assert list(diff) == [c], f"{m} {st}: alleles differ at {diff}, expected only {c}"

        w = np.arange(c - a.contrib_flank, c + a.contrib_flank + 1)
        for ai, alab in enumerate(["REF", "ALT"]):
            for j in w:
                onehot = seq[ai, :, j]
                for k in range(4):
                    contrib_rows.append((m, st, alab, j - c, NUC[k], int(onehot[k] == 1),
                                         hyp[ai, k, j], hyp[ai, k, j] * onehot[k]))

    # ---------------- predictions ----------------
    tracks = {"REF": [], "ALT": []}
    for f in folds:
        pre = os.path.join(pf, f"{m}_fold{f}")
        with h5py.File(pre + ".variant_predictions.h5", "r") as h:
            o = h["observed"]
            for alab, key in [("REF", "allele1"), ("ALT", "allele2")]:
                logits = np.asarray(o[f"{key}_pred_profiles"])[0].ravel()   # fwd/rc-averaged logits
                counts = float(np.asarray(o[f"{key}_pred_counts"]).ravel()[0])
                prob = np.exp(logits - logits.max()); prob /= prob.sum()
                tracks[alab].append(prob * counts)
        s = pd.read_csv(pre + ".variant_scores.tsv", sep="\t").iloc[0]
        score_rows.append((m, f, s.get("logfc"), s.get("jsd"),
                           s.get("allele1_pred_counts"), s.get("allele2_pred_counts")))
    for alab in ["REF", "ALT"]:
        tr = np.mean(tracks[alab], axis=0)
        out_len = len(tr); oc = out_len // 2                 # SNP sits at output index out_len//2
        fl = min(a.track_flank, oc - 1)
        for j in range(oc - fl, oc + fl + 1):
            pred_rows.append((m, alab, j - oc, tr[j]))

# ---------------- bigwig tracks ----------------
bw_rows = []
for m in models:
    base = f"{a.cb}/ATAC_{m}_merged_alldonors/biasFactor{a.bias}"
    # observed signal: pipeline's Tn5-shifted training bigWig, same source for all models
    cands = sorted(glob.glob(f"{base}/fold0/bias_model/auxiliary/*_data_unstranded.bw"))
    obs = cands[0] if cands else f"{base}/data/ATAC_{m}_merged_alldonors.actual.bw"
    print(f"{m} observed bigWig: {obs}")
    bws = {"observed": obs}
    if a.pred_bw_pattern:
        bws["predicted_bw"] = a.pred_bw_pattern.format(cb=a.cb, model=m, bias=a.bias)
    for tname, fn in bws.items():
        if not os.path.exists(fn):
            print(f"WARNING: missing {fn}, skipped"); continue
        bw = pyBigWig.open(fn)
        start, end = snp0 - a.track_flank, snp0 + a.track_flank + 1
        vals = np.nan_to_num(np.array(bw.values(chrom, start, end), dtype=float))
        # depth normalisation: genome-wide sum of the track (total Tn5 cuts) -> cuts per million,
        # so tracks from libraries of different depth can share one y-axis
        total = float(bw.header().get("sumData", 0.0))
        bw.close()
        scale = 1e6 / total if total > 0 else np.nan
        print(f"{m} {tname}: genome-wide sum {total:.4g} -> CPM scale {scale:.4g}")
        for i, v in enumerate(vals):
            bw_rows.append((m, tname, start + i - snp0, v, v * scale))

# ---------------- ISM (optional: per_fold/<model>.ism.tsv) ----------------
ism_tabs, ism_summ = [], []
for m in models:
    fn = os.path.join(pf, f"{m}.ism.tsv")
    if not os.path.exists(fn):
        print(f"WARNING: no ISM for {m} ({fn}), skipped"); continue
    d = pd.read_csv(fn, sep="\t")
    g = (d.groupby(["pos_rel", "ref_base", "mut_base"], as_index=False)
          .agg(mean_log2fc=("log2fc", "mean"), sd_log2fc=("log2fc", "std")))
    g.insert(0, "model", m)
    g["is_variant"] = ((g.pos_rel == 0) & (g.mut_base == alt)).astype(int)
    # rank of each substitution by |effect| (1 = most disruptive in the window)
    g["rank_abs"] = g["mean_log2fc"].abs().rank(ascending=False, method="min").astype(int)
    # rank of each position by its strongest substitution
    pos_max = g.groupby("pos_rel")["mean_log2fc"].apply(lambda x: x.abs().max())
    pos_rank = pos_max.rank(ascending=False, method="min").astype(int)
    g["pos_rank"] = g["pos_rel"].map(pos_rank)
    ism_tabs.append(g)
    vr = g[g.is_variant == 1]
    if len(vr) == 1:
        ism_summ.append(dict(model=m, ism_log2fc=float(vr.mean_log2fc.iloc[0]),
                             ism_rank=int(vr.rank_abs.iloc[0]), ism_n=len(g),
                             ism_pos_rank=int(vr.pos_rank.iloc[0]), ism_n_pos=int(g.pos_rel.nunique())))

# ---------------- write ----------------
o = lambda s: os.path.join(a.outdir, f"{a.rsid}.{s}")
pd.DataFrame(contrib_rows, columns=["model", "shap_type", "allele", "pos_rel", "nucleotide",
                                    "is_ref_base", "hyp", "proj"]).to_csv(o("contrib_window.tsv"), sep="\t", index=False)
pd.DataFrame(pred_rows, columns=["model", "allele", "pos_rel", "pred"]).to_csv(o("pred_profiles.tsv"), sep="\t", index=False)
pd.DataFrame(bw_rows, columns=["model", "track", "pos_rel", "value", "value_cpm"]).to_csv(o("bigwig_tracks.tsv"), sep="\t", index=False)
scores = pd.DataFrame(score_rows, columns=["model", "fold", "logfc", "jsd",
                                           "allele1_pred_counts", "allele2_pred_counts"])
scores.to_csv(o("scores_by_fold.tsv"), sep="\t", index=False)

# per-model summary across folds (logfc = log2(ALT/REF))
summ = (scores.groupby("model", sort=False)
        .agg(mean_log2fc=("logfc", "mean"), sd_log2fc=("logfc", "std"),
             folds_pos=("logfc", lambda x: int((x > 0).sum())), n_folds=("logfc", "size"),
             mean_jsd=("jsd", "mean"))
        .reset_index())
summ.insert(3, "alt_over_ref", 2 ** summ["mean_log2fc"])
if ism_summ:
    summ = summ.merge(pd.DataFrame(ism_summ), on="model", how="left")
    pd.concat(ism_tabs).to_csv(o("ism.tsv"), sep="\t", index=False, float_format="%.5g")
summ.to_csv(o("scores_summary.tsv"), sep="\t", index=False, float_format="%.4g")
print(summ.to_string(index=False))
print("Wrote", a.outdir)

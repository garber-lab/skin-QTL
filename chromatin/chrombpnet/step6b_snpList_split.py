#!/usr/bin/env python
"""
Split batched (many-variant) ChromBPNet outputs into the per-variant layout that the step6
scripts expect, so step6_summarize / step6_motif_annotate / step6_plot run unchanged.

Batched inputs, per model (written by step6b_snpList_run.sh):
  <batch>/per_fold/<model>_fold<f>.variant_shap.{counts,profile}.h5   (2N regions: REF rows then ALT)
  <batch>/per_fold/<model>_fold<f>.variant_predictions.h5             (N rows, input order)
  <batch>/per_fold/<model>_fold<f>.variant_scores.tsv                 (N rows, variant_id column)
  <batch>/per_fold/<model>.ism.tsv                                    (variant_id column)
  <batch>/motifs/finemo/<model>/<head>/hits.tsv                       (peak_id 0..2N-1)

Per-variant outputs, under <outroot>/<rsID>/:
  <rsID>.variant.tsv
  per_fold/<model>_fold<f>.variant_shap.{counts,profile}.h5
  per_fold/<model>_fold<f>.variant_predictions.h5
  per_fold/<model>_fold<f>.variant_scores.tsv
  per_fold/<model>.ism.tsv
  motifs/finemo/<model>/<head>/hits.tsv  (peak_id rewritten to 0 = REF, 1 = ALT) + modisco_source.txt
"""
import argparse, os, shutil, sys
import numpy as np
import pandas as pd
import h5py
import deepdish as dd

p = argparse.ArgumentParser()
p.add_argument("--batch", required=True, help="batch working directory")
p.add_argument("--outroot", required=True, help="posterchild directory")
p.add_argument("--variants", required=True, help="chrombpnet-schema table, no header")
p.add_argument("--models", required=True)
p.add_argument("--folds", default="0,1,2,3,4")
p.add_argument("--heads", default="counts,profile")
a = p.parse_args()

models, folds, heads = a.models.split(","), a.folds.split(","), a.heads.split(",")
var = pd.read_csv(a.variants, sep="\t", header=None, names=["chr", "pos", "ref", "alt", "id"])
ids = var["id"].astype(str).tolist()
N = len(ids)
print(f"splitting {N} variants x {len(models)} models x {len(folds)} folds", flush=True)

# per-variant folders + variant.tsv
for i, rsid in enumerate(ids):
    d = os.path.join(a.outroot, rsid)
    os.makedirs(os.path.join(d, "per_fold"), exist_ok=True)
    os.makedirs(os.path.join(d, "logs"), exist_ok=True)
    var.iloc[[i]].to_csv(os.path.join(d, f"{rsid}.variant.tsv"), sep="\t", header=False, index=False)

pf = os.path.join(a.batch, "per_fold")
missing = []

for m in models:
    for f in folds:
        # ---------- contribution scores: rows i (REF) and N+i (ALT) ----------
        for h in heads:
            fn = os.path.join(pf, f"{m}_fold{f}.variant_shap.{h}.h5")
            if not os.path.exists(fn):
                missing.append(fn); continue
            d = dd.io.load(fn)
            vids = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in np.asarray(d["variant_ids"])])
            alleles = np.asarray(d["alleles"])
            for i, rsid in enumerate(ids):
                sel = np.where(vids == rsid)[0]
                if len(sel) != 2:
                    missing.append(f"{fn}:{rsid}"); continue
                sel = sel[np.argsort(alleles[sel], kind="stable")]     # REF row first, then ALT
                out = {k: {"seq": np.asarray(d[k]["seq"])[sel]} for k in ("raw", "shap", "projected_shap") if k in d}
                out["variant_ids"] = vids[sel]
                out["alleles"] = alleles[sel]
                dd.io.save(os.path.join(a.outroot, rsid, "per_fold",
                                        f"{m}_fold{f}.variant_shap.{h}.h5"), out)

        # ---------- predictions: row i ----------
        fn = os.path.join(pf, f"{m}_fold{f}.variant_predictions.h5")
        if os.path.exists(fn):
            with h5py.File(fn, "r") as hf:
                obs = {k: np.asarray(hf["observed"][k]) for k in hf["observed"].keys()}
            for i, rsid in enumerate(ids):
                o = os.path.join(a.outroot, rsid, "per_fold", f"{m}_fold{f}.variant_predictions.h5")
                with h5py.File(o, "w") as ho:
                    g = ho.create_group("observed")
                    for k, arr in obs.items():
                        g.create_dataset(k, data=arr[i:i + 1])
        else:
            missing.append(fn)

        # ---------- scores table: row for this variant ----------
        fn = os.path.join(pf, f"{m}_fold{f}.variant_scores.tsv")
        if os.path.exists(fn):
            sc = pd.read_csv(fn, sep="\t")
            key = "variant_id" if "variant_id" in sc.columns else sc.columns[4]
            for rsid in ids:
                sub = sc[sc[key].astype(str) == rsid]
                if len(sub) != 1:
                    missing.append(f"{fn}:{rsid}"); continue
                sub.to_csv(os.path.join(a.outroot, rsid, "per_fold",
                                        f"{m}_fold{f}.variant_scores.tsv"), sep="\t", index=False)
        else:
            missing.append(fn)

    # ---------- ISM ----------
    fn = os.path.join(pf, f"{m}.ism.tsv")
    if os.path.exists(fn):
        ism = pd.read_csv(fn, sep="\t")
        for rsid, g in ism.groupby("variant_id"):
            if rsid not in ids:
                continue
            g.drop(columns=["variant_id"]).to_csv(
                os.path.join(a.outroot, str(rsid), "per_fold", f"{m}.ism.tsv"), sep="\t", index=False)
    else:
        missing.append(fn)

    # ---------- Fi-NeMo hits: peak_id 0..N-1 = REF, N..2N-1 = ALT ----------
    for h in heads:
        hd = os.path.join(a.batch, "motifs", "finemo", m, h)
        fn = os.path.join(hd, "hits.tsv")
        if not os.path.exists(fn):
            missing.append(fn); continue
        hits = pd.read_csv(fn, sep="\t")
        hits["vi"] = hits["peak_id"] % N
        hits["allele_idx"] = (hits["peak_id"] // N).astype(int)     # 0 = REF, 1 = ALT
        for vi, g in hits.groupby("vi"):
            rsid = ids[int(vi)]
            od = os.path.join(a.outroot, rsid, "motifs", "finemo", m, h)
            os.makedirs(od, exist_ok=True)
            g = g.drop(columns=["vi"]).copy()
            g["peak_id"] = g.pop("allele_idx")                      # per-variant convention
            g.to_csv(os.path.join(od, "hits.tsv"), sep="\t", index=False)
            src = os.path.join(hd, "modisco_source.txt")
            if os.path.exists(src):
                shutil.copyfile(src, os.path.join(od, "modisco_source.txt"))
        # variants with no hits still need the folder (and an empty table) for the annotation step
        for vi, rsid in enumerate(ids):
            od = os.path.join(a.outroot, rsid, "motifs", "finemo", m, h)
            if not os.path.exists(os.path.join(od, "hits.tsv")):
                os.makedirs(od, exist_ok=True)
                hits.iloc[0:0].drop(columns=["vi", "allele_idx"]).to_csv(
                    os.path.join(od, "hits.tsv"), sep="\t", index=False)
            src = os.path.join(hd, "modisco_source.txt")
            if os.path.exists(src) and not os.path.exists(os.path.join(od, "modisco_source.txt")):
                shutil.copyfile(src, os.path.join(od, "modisco_source.txt"))
    print(f"  {m} done", flush=True)

if missing:
    print(f"WARNING: {len(missing)} missing/short inputs, first 10:", file=sys.stderr)
    for x in missing[:10]:
        print("   ", x, file=sys.stderr)
print(f"Wrote per-variant outputs under {a.outroot}")

#!/usr/bin/env python
"""
Write a Fi-NeMo name file (-N) for EVERY MoDISco pattern of each model, so hit tables carry
readable motif names instead of pos_patterns.pattern_7.

Pattern numbering is per model, so each model gets its own file; the same TF may be pattern_2 in
one model and pattern_6 in another. Names come from the TOMTOM matches in that model's motifs.html.

Usage:
  python step7_motifScan_names.py --cb /path/chrombpnet --models FRB_IFN,... --head profile \
      --out /path/motifscan/names
Outputs per model, under <out>/<model>/:
  names.txt        pattern <tab> <TF>_<pattern id>      (for finemo call-hits -N)
  patterns.tsv     pattern, num_seqlets, top match, qval, cleaned TF name, sign
"""
import argparse, os, re, sys
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--cb", required=True)
p.add_argument("--bias", default="0.5")
p.add_argument("--models", required=True)
p.add_argument("--head", default="profile", choices=["profile", "counts"])
p.add_argument("--min_seqlets", type=int, default=0, help="0 keeps every pattern")
p.add_argument("--out", required=True)
a = p.parse_args()

REPORT = {"profile": "modisco_profile", "counts": "modisco_counts"}[a.head]


def parse_motifs_html(fn):
    x = open(fn, encoding="utf-8", errors="replace").read()
    rows = re.split(r"<tr", x)[1:]
    cells = [[re.sub(r"<[^>]+>", "", c).strip()
              for c in re.findall(r"(?s)<t[hd][^>]*>(.*?)</t[hd]>", r)] for r in rows]
    if not cells:
        return None
    hdr, body = cells[0], [c for c in cells[1:] if len(c) == len(cells[0])]
    return pd.DataFrame(body, columns=hdr) if body else None


def clean_tf(x):
    """FOSL1+JUNB_MA1137.1 -> FOSL1+JUNB ; Tp73.mouse_p53l_1 -> TP73 ; ATF2_HUMAN.H11MO.0.B -> ATF2"""
    if not x:
        return "NA"
    y = re.sub(r"_.*$", "", x)
    y = re.sub(r"\.(mouse|human)$", "", y, flags=re.I)
    return y.upper() or "NA"


all_rows = []
for m in a.models.split(","):
    rep = os.path.join(a.cb, f"ATAC_{m}_merged_alldonors", f"biasFactor{a.bias}", "fold0",
                       "chrombpnet_model", "evaluation", REPORT, "motifs.html")
    if not os.path.exists(rep):
        print(f"WARNING: no report for {m}: {rep}", file=sys.stderr); continue
    df = parse_motifs_html(rep)
    if df is None or "pattern" not in df.columns:
        print(f"WARNING: could not parse {rep}", file=sys.stderr); continue

    d = os.path.join(a.out, m)
    os.makedirs(d, exist_ok=True)
    rows, used = [], {}
    with open(os.path.join(d, "names.txt"), "w") as nf:
        for _, r in df.iterrows():
            pat = str(r["pattern"])
            try:
                seqlets = int(float(r.get("num_seqlets", 0) or 0))
            except ValueError:
                seqlets = 0
            if seqlets < a.min_seqlets:
                continue
            tf = clean_tf(str(r.get("match0", "")))
            pid = pat.split(".")[-1].replace("pattern_", "p")
            sign = "neg" if pat.startswith("neg") else "pos"
            name = f"{tf}_{pid}" if sign == "pos" else f"{tf}_{pid}_neg"
            used[name] = used.get(name, 0) + 1
            if used[name] > 1:                      # keep names unique within a model
                name = f"{name}_{used[name]}"
            nf.write(f"{pat}\t{name}\n")
            rows.append(dict(model=m, pattern=pat, name=name, tf=tf, sign=sign,
                             num_seqlets=seqlets, match0=r.get("match0", ""), qval0=r.get("qval0", "")))
    pd.DataFrame(rows).to_csv(os.path.join(d, "patterns.tsv"), sep="\t", index=False)
    all_rows += rows
    print(f"{m}: {len(rows)} patterns -> {d}/names.txt")

if not all_rows:
    sys.exit("ERROR: no patterns found in any model")
pd.DataFrame(all_rows).to_csv(os.path.join(a.out, f"patterns_all_models_{a.head}.tsv"), sep="\t", index=False)
print(f"Wrote {os.path.join(a.out, f'patterns_all_models_{a.head}.tsv')} ({len(all_rows)} patterns)")

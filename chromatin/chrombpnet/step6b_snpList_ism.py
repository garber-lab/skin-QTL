#!/usr/bin/env python
"""
In-silico saturation mutagenesis for MANY variants, all folds of one model.
Batched version of step6_ism_snp.py: each fold's model is loaded once and every variant scored
against it, instead of reloading the model per variant.

For each variant, every position within +/- flank bp is substituted to each of the 3 non-reference
bases in the REFERENCE sequence; log2FC = log2(mutant counts / REF counts), counts averaged over
forward and reverse-complement strands (same convention as variant_scoring.py).

Run inside the chrombpnet container:
  python step6b_snpList_ism.py --vs_src /path/variant-scorer/src \
      --variants variants.tsv --genome hg38.fa \
      --model_glob '/path/.../fold{f}/chrombpnet_model/models/*chrombpnet_nobias.h5' \
      --folds 0,1,2,3,4 --flank 50 --out batch/KRT_IFN.ism.tsv

variants.tsv: chr, pos(1-based), REF, ALT, variant_id - no header (chrombpnet schema)
Output TSV: variant_id, fold, pos_rel, ref_base, mut_base, log2fc
"""
import argparse, glob, os, sys, time

p = argparse.ArgumentParser()
p.add_argument("--vs_src", default=os.getcwd(), help="variant-scorer/src (for its model loader)")
p.add_argument("--variants", required=True)
p.add_argument("--genome", required=True)
p.add_argument("--model_glob", required=True, help="path pattern with {f} for the fold number")
p.add_argument("--folds", default="0,1,2,3,4")
p.add_argument("--flank", type=int, default=50)
p.add_argument("--bs", type=int, default=256)
p.add_argument("--out", required=True)
a = p.parse_args()

sys.path.insert(0, a.vs_src)          # python puts the script's folder on sys.path, not the cwd
import numpy as np
import pandas as pd
import pyfaidx
import tensorflow as tf
from utils.helpers import load_model_wrapper
from utils.one_hot import dna_to_one_hot

var = pd.read_csv(a.variants, sep="\t", header=None, names=["chr", "pos", "ref", "alt", "id"])
genome = pyfaidx.Fasta(a.genome)
BASES = "ACGT"
print(f"{len(var)} variants, folds {a.folds}, flank +/-{a.flank}", flush=True)

rows = []
for f in a.folds.split(","):
    hits = sorted(glob.glob(a.model_glob.format(f=f)))
    if len(hits) != 1:
        sys.exit(f"ERROR: fold {f}: expected 1 model for {a.model_glob.format(f=f)}, found {len(hits)}")
    model = load_model_wrapper(hits[0])
    L = model.input_shape[1] if not isinstance(model.input_shape, list) else model.input_shape[0][1]
    half = L // 2
    t0 = time.time()

    for _, v in var.iterrows():                      # one model load, every variant scored against it
        chrom, pos0, ref = str(v.chr), int(v.pos) - 1, str(v.ref).upper()
        seq = str(genome[chrom][pos0 - half:pos0 + half]).upper()
        if len(seq) != L:
            print(f"WARNING: {v.id}: sequence length {len(seq)} != {L}, skipped", file=sys.stderr); continue
        if seq[half] != ref:
            # VCF REF/ALT swapped relative to hg38: substitute, as variant_shap/variant_scoring do
            print(f"WARNING: {v.id}: genome base {seq[half]} != VCF REF {ref}; substituting REF", file=sys.stderr)
            seq = seq[:half] + ref + seq[half + 1:]

        seqs, meta = [seq], [(None, None, None)]
        for j in range(-a.flank, a.flank + 1):
            i = half + j
            rb = seq[i]
            for b in BASES:
                if b == rb:
                    continue
                seqs.append(seq[:i] + b + seq[i + 1:])
                meta.append((j, rb, b))

        x = dna_to_one_hot(seqs)
        fwd = model.predict(x, batch_size=a.bs, verbose=0)[1].ravel()
        rc = model.predict(x[:, ::-1, ::-1], batch_size=a.bs, verbose=0)[1].ravel()
        counts = (np.exp(fwd) + np.exp(rc)) / 2
        l2 = np.log2(counts / counts[0])
        for (j, rb, b), val in zip(meta[1:], l2[1:]):
            rows.append((v.id, f, j, rb, b, float(val)))

    print(f"fold {f}: {len(var)} variants in {time.time() - t0:.0f}s ({hits[0]})", flush=True)
    del model
    tf.keras.backend.clear_session()

os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
pd.DataFrame(rows, columns=["variant_id", "fold", "pos_rel", "ref_base", "mut_base", "log2fc"]) \
  .to_csv(a.out, sep="\t", index=False)
print("Wrote", a.out)

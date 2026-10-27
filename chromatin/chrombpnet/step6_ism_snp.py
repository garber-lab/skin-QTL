#!/usr/bin/env python
"""
In-silico saturation mutagenesis (ISM) around one variant, for all folds of one model.

For every position within +/- flank bp of the variant, each of the 3 non-reference
bases is substituted into the REFERENCE sequence and the model predicts total counts.
log2FC = log2(mutant counts / REF counts), with counts averaged over forward and
reverse-complement strands (same convention as variant-scorer's variant_scoring.py,
so the G->A value at pos_rel 0 should reproduce its logfc).

Run inside the container:
  python step6_ism_snp.py --vs_src /path/variant-scorer/src --variant rsX.variant.tsv --genome hg38.fa \
      --model_glob '/path/ATAC_KRT_IFN_merged_alldonors/biasFactor0.5/fold{f}/chrombpnet_model/models/*chrombpnet_nobias.h5' \
      --folds 0,1,2,3,4 --flank 50 --out per_fold/KRT_IFN.ism.tsv

Output TSV: fold, pos_rel, ref_base, mut_base, log2fc
"""
import argparse, glob, os, sys

p = argparse.ArgumentParser()
p.add_argument("--vs_src", default=os.getcwd(), help="variant-scorer/src (for its model loader)")
p.add_argument("--variant", required=True)
p.add_argument("--genome", required=True)
p.add_argument("--model_glob", required=True, help="path pattern with {f} for the fold number")
p.add_argument("--folds", default="0,1,2,3,4")
p.add_argument("--flank", type=int, default=50)
p.add_argument("--bs", type=int, default=128)
p.add_argument("--out", required=True)
a = p.parse_args()

# python puts the script's own folder on sys.path, not the cwd, so add variant-scorer/src explicitly
sys.path.insert(0, a.vs_src)
import numpy as np
import pandas as pd
import pyfaidx
import tensorflow as tf
from utils.helpers import load_model_wrapper
from utils.one_hot import dna_to_one_hot

v = pd.read_csv(a.variant, sep="\t", header=None, names=["chr", "pos", "ref", "alt", "id"]).iloc[0]
chrom, pos0, ref = str(v.chr), int(v.pos) - 1, v.ref.upper()
genome = pyfaidx.Fasta(a.genome)
BASES = "ACGT"

rows = []
for f in a.folds.split(","):
    hits = sorted(glob.glob(a.model_glob.format(f=f)))
    if len(hits) != 1:
        sys.exit(f"ERROR: fold {f}: expected 1 model for {a.model_glob.format(f=f)}, found {len(hits)}")
    model = load_model_wrapper(hits[0])
    L = model.input_shape[1] if not isinstance(model.input_shape, list) else model.input_shape[0][1]
    half = L // 2

    # same centering as variant-scorer's VariantGenerator: SNP at index L//2
    seq = str(genome[chrom][pos0 - half:pos0 + half]).upper()
    assert len(seq) == L, f"sequence length {len(seq)} != {L}"
    if seq[half] != ref:
        # VCF REF/ALT can be swapped relative to hg38; variant-scorer substitutes allele1 regardless,
        # so do the same here and keep the REF/ALT labels as the VCF gives them
        print(f"WARNING: genome base {seq[half]} != VCF REF {ref}; substituting REF", file=sys.stderr)
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
        rows.append((f, j, rb, b, float(val)))
    print(f"fold {f}: {len(seqs) - 1} substitutions scored ({hits[0]})")

    del model
    tf.keras.backend.clear_session()

pd.DataFrame(rows, columns=["fold", "pos_rel", "ref_base", "mut_base", "log2fc"]).to_csv(a.out, sep="\t", index=False)
print("Wrote", a.out)

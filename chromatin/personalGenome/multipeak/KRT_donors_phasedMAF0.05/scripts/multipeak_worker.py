"""
Personal chromBPNet prediction worker (runs inside the chrombpnet container).
Called by multipeak_chunk.sh in two modes:

  predict   <chunk.tsv> <ref_windows.fa> <genotypes.tsv> <samples.txt> <chunk_dir> <input_len> <batch_size>
            <strict_ref> <id_map|none> <n_threads> <group> LABEL=model.h5 [LABEL=model.h5 ...]
      Build every donor's two haplotypes per peak (biallelic SNVs), collect the unique haplotypes of all
      peaks in the chunk, predict them with this job's models, save
        <chunk_dir>/structure_g<group>.pkl  (haplotype structure; identical for every model group)
        <chunk_dir>/pred/<LABEL>.npy        (predicted counts per unique haplotype)
        <chunk_dir>/timing_g<group>.tsv

  finalize  <chunk_dir> <save_folds> <conditions> <labels>
      Once all model groups of the chunk are done: fold-averaged columns per condition, then write
        peak_info.tsv, haplotypes.tsv.gz, donor_hap_ids.tsv.gz, diploid_<column>.tsv.gz
"""
import gzip
import os
import pickle
import sys
import time

import numpy as np

LUT = {b: i for i, b in enumerate("ACGT")}


def read_haplotypes(bed, ref_fa, geno_tsv, samples_txt, input_len, strict_ref, id_map):
    import bisect

    # --- peaks and REF sequences (same order) ---
    peaks = []
    for line in open(bed):
        f = line.rstrip("\n").split("\t")
        peaks.append(dict(chrom=f[0], start=int(f[1]), end=int(f[2]), name=f[3], summit=int(f[9]),
                          center=int(f[10]), ws=int(f[11])))
    seqs, cur = [], None
    for line in open(ref_fa):
        if line.startswith(">"):
            if cur is not None:
                seqs.append("".join(cur).upper())
            cur = []
        else:
            cur.append(line.strip())
    seqs.append("".join(cur).upper())
    if len(seqs) != len(peaks):
        sys.exit(f"ERROR {len(seqs)} REF sequences for {len(peaks)} peaks")
    for p, s in zip(peaks, seqs):
        if len(s) != input_len:
            sys.exit(f"ERROR {p['name']}: REF length {len(s)} != {input_len}")
        p["seq"] = s

    samples = [l.strip() for l in open(samples_txt) if l.strip()]
    if id_map != "none":                        # rename VCF sample IDs to donor names for the output headers
        rename = dict(l.split()[:2] for l in open(id_map) if len(l.split()) >= 2)
        missing_ids = [s for s in samples if s not in rename]
        if missing_ids:
            sys.exit(f"ERROR {len(missing_ids)} VCF samples not in ID_MAP, e.g. {missing_ids[:3]}")
        samples = [rename[s] for s in samples]
    n_samples = len(samples)

    # --- variant lines grouped by chromosome, sorted by position, de-duplicated ---
    by_chrom = {}
    for line in open(geno_tsv):
        chrom, pos, vid, ref, alt, rest = line.split("\t", 5)
        by_chrom.setdefault(chrom, {})[(int(pos), ref, alt)] = line     # overlapping windows -> one copy
    for c in by_chrom:
        items = sorted(by_chrom[c].items())
        by_chrom[c] = ([k[0] for k, _ in items], [l for _, l in items])

    # --- per peak: donor haplotypes, unique haplotype list ---
    n_entries = 0
    for p in peaks:
        pos_list, lines = by_chrom.get(p["chrom"], ([], []))
        lo = bisect.bisect_left(pos_list, p["ws"] + 1)
        hi = bisect.bisect_left(pos_list, p["ws"] + input_len + 1)
        variants, hap_vars = [], [[[], []] for _ in range(n_samples)]
        n_indel = n_multi = n_refmis = 0
        n_het = np.zeros(n_samples, int)             # heterozygous SNVs per donor in this window
        n_het_unphased = np.zeros(n_samples, int)    # ... of which written without phase ("/")
        missing = np.zeros(n_samples, bool)
        for line in lines[lo:hi]:
            f = line.rstrip("\n").split("\t")
            pos, ref, alt, vtype, gts = int(f[1]), f[3].upper(), f[4].upper(), f[5], f[6:]
            if "," in alt:
                n_multi += 1
                continue
            if vtype != "SNP" or len(ref) != 1 or alt not in LUT:
                n_indel += 1
                continue
            if len(gts) != n_samples:
                sys.exit(f"ERROR {len(gts)} genotypes at {p['chrom']}:{pos}, expected {n_samples}")
            off = pos - 1 - p["ws"]
            if p["seq"][off] != ref:
                if strict_ref == "1":
                    sys.exit(f"ERROR REF mismatch at {p['chrom']}:{pos}: VCF {ref} vs FASTA {p['seq'][off]} "
                             "(wrong genome build or strand? set STRICT_REF=0 to skip such sites)")
                n_refmis += 1
                continue
            v = len(variants)
            variants.append((pos, ref, alt, off))
            for s, gt in enumerate(gts):
                a = gt.replace("/", "|").split("|")
                if len(a) == 1:
                    a = a * 2
                if "." in a:
                    missing[s] = True                # missing allele -> left as REF
                al = [x == "1" for x in a]
                if al[0] != al[1]:
                    n_het[s] += 1
                    if "|" not in gt:
                        n_het_unphased[s] += 1       # het without phase: order taken as written
                for h in (0, 1):
                    if al[h]:
                        hap_vars[s][h].append(v)
        hap_index = {(): 0}                          # local hap_id 0 = hg38 REF
        hap_ids = np.zeros((n_samples, 2), dtype=np.int32)
        for s in range(n_samples):
            for h in (0, 1):
                key = tuple(hap_vars[s][h])
                if key not in hap_index:
                    hap_index[key] = len(hap_index)
                hap_ids[s, h] = hap_index[key]
        keys = sorted(hap_index, key=hap_index.get)
        p.update(variants=variants, hap_ids=hap_ids, keys=keys, offset=n_entries,
                 n_indel=n_indel, n_multi=n_multi, n_refmis=n_refmis, n_missing=int(missing.sum()),
                 # phase only changes the diploid prediction if a donor has >= 2 hets and not all are phased
                 n_ambiguous=int(((n_het >= 2) & (n_het_unphased >= 1)).sum()),
                 n_carriers=np.bincount(hap_ids.ravel(), minlength=len(keys)))
        n_entries += len(keys)
    return peaks, samples, n_entries


def encode(p, key, ref_oh):
    x = ref_oh.copy()
    for v in key:
        off, alt = p["variants"][v][3], p["variants"][v][2]
        x[off] = 0.0
        x[off, LUT[alt]] = 1.0
    return x


def predict(args):
    (bed, ref_fa, geno_tsv, samples_txt, chunk_dir, input_len, batch_size,
     strict_ref, id_map, n_threads, group) = args[:11]
    models = [m.split("=", 1) for m in args[11:]]          # [(label, path)], label = CONDITION_foldN
    input_len, batch_size, n_threads = int(input_len), int(batch_size), int(n_threads)

    t0 = time.time()
    import tensorflow as tf
    # use exactly the cores LSF gave this job (TF otherwise sizes its thread pools to the whole node)
    tf.config.threading.set_intra_op_parallelism_threads(n_threads)
    tf.config.threading.set_inter_op_parallelism_threads(min(2, n_threads))
    timing = [("tf_import", time.time() - t0, 0)]
    print(f"TensorFlow {tf.__version__}, numpy {np.__version__}, threads {n_threads}, "
          f"GPUs: {len(tf.config.list_physical_devices('GPU'))}")

    t1 = time.time()
    peaks, samples, n_entries = read_haplotypes(bed, ref_fa, geno_tsv, samples_txt, input_len, strict_ref, id_map)
    timing.append(("parse_haplotypes", time.time() - t1, n_entries))
    print(f"Peaks: {len(peaks)}; samples: {len(samples)}; unique haplotypes (incl. REF per peak): {n_entries}")

    # one-hot REF per peak (A=1000 C=0100 G=0010 T=0001 N=0000) and the global haplotype list
    ref_oh = []
    entries = []
    for pi, p in enumerate(peaks):
        oh = np.zeros((input_len, 4), dtype=np.float32)
        for i, b in enumerate(p["seq"]):
            if b in LUT:
                oh[i, LUT[b]] = 1.0
        ref_oh.append(oh)
        entries += [(pi, k) for k in p["keys"]]

    # haplotype structure for the finalize step. Identical for every model group, but each group writes
    # its OWN file: groups of a chunk start together, and replacing one shared file at the same time
    # can fail on a network filesystem.
    struct_path = os.path.join(chunk_dir, f"structure_g{group}.pkl")
    tmp = f"{struct_path}.tmp"
    keep = ("chrom", "start", "end", "name", "summit", "center", "ws", "variants", "hap_ids", "keys",
            "offset", "n_indel", "n_multi", "n_refmis", "n_missing", "n_ambiguous", "n_carriers")
    with open(tmp, "wb") as fh:
        pickle.dump(dict(samples=samples, input_len=input_len, n_entries=n_entries,
                         peaks=[{k: p[k] for k in keep} for p in peaks]), fh)
    os.replace(tmp, struct_path)

    os.makedirs(os.path.join(chunk_dir, "pred"), exist_ok=True)
    for lab, path in models:
        t = time.time()
        try:
            model = tf.keras.models.load_model(path, compile=False)
        except Exception as e:
            print(f"Plain load failed ({e}); retrying with custom_objects")
            model = tf.keras.models.load_model(path, compile=False,
                                               custom_objects={"tf": tf, "multinomial_nll": lambda *a: 0.0})
        if model.input_shape[1] != input_len:
            sys.exit(f"ERROR {lab}: model input length {model.input_shape[1]} != {input_len}")
        timing.append((f"load_{lab}", time.time() - t, 0))

        # fixed batch size (last batch padded with zeros) -> the model's graph is built only once
        t = time.time()
        lc_all = []
        for b0 in range(0, n_entries, batch_size):
            batch = entries[b0:b0 + batch_size]
            x = np.zeros((batch_size, input_len, 4), dtype=np.float32)
            for j, (pi, key) in enumerate(batch):
                x[j] = encode(peaks[pi], key, ref_oh[pi])
            out = model.predict_on_batch(x)
            lc_all.append(np.asarray(out[1]).reshape(-1)[:len(batch)])
        # ChromBPNet's count head is trained on log(1 + total reads), so expm1 converts it back to reads
        counts = np.expm1(np.concatenate(lc_all)) if lc_all else np.zeros(0)
        tmp = os.path.join(chunk_dir, "pred", f"{lab}.tmp.npy")
        np.save(tmp, counts)
        os.replace(tmp, os.path.join(chunk_dir, "pred", f"{lab}.npy"))
        dt = time.time() - t
        timing.append((f"predict_{lab}", dt, n_entries))
        print(f"  {lab}: load {timing[-2][1]:.1f} s, predict {dt:.1f} s ({n_entries / max(dt, 1e-9):.1f} haplotypes/s)")
        tf.keras.backend.clear_session()
        del model

    chunk_id = os.path.basename(os.path.normpath(chunk_dir)).replace("chunk_", "")
    with open(os.path.join(chunk_dir, f"timing_g{group}.tsv"), "w") as fh:
        for step, sec, n in timing:
            fh.write(f"{chunk_id}\t{group}\t{step}\t{sec:.2f}\t{n}\n")


def finalize(args):
    chunk_dir, save_folds, conditions, labels = args[:4]
    conditions, labels = conditions.split(), labels.split()
    t = time.time()
    import glob
    structs = sorted(glob.glob(os.path.join(chunk_dir, "structure_g*.pkl")))
    if not structs:
        sys.exit(f"ERROR no structure_g*.pkl in {chunk_dir}")
    st = pickle.load(open(structs[0], "rb"))         # all groups wrote the same structure; any one will do
    peaks, samples, input_len = st["peaks"], st["samples"], st["input_len"]

    counts = {lab: np.load(os.path.join(chunk_dir, "pred", f"{lab}.npy")) for lab in labels}
    for lab in labels:
        if len(counts[lab]) != st["n_entries"]:
            sys.exit(f"ERROR {lab}: {len(counts[lab])} predictions, expected {st['n_entries']}")
    # fold-averaged counts per condition (mean of the folds' counts, per haplotype)
    for c in conditions:
        counts[f"{c}_mean"] = np.mean([counts[l] for l in labels if l.rsplit("_fold", 1)[0] == c], axis=0)
    cols = []
    for c in conditions:
        cols += [l for l in labels if l.rsplit("_fold", 1)[0] == c] + [f"{c}_mean"]
    matrix_cols = cols if save_folds == "1" else [f"{c}_mean" for c in conditions]

    with open(os.path.join(chunk_dir, "peak_info.tsv"), "w") as fh:
        fh.write("peak_id\tchrom\tstart\tend\tsummit_offset\tcenter\twindow\tn_snvs\tn_skipped_indel\t"
                 "n_skipped_multiallelic\tn_skipped_ref_mismatch\tn_unique_haps\tn_donors_phase_ambiguous\t"
                 "n_donors_missing\n")
        for p in peaks:
            fh.write(f"{p['name']}\t{p['chrom']}\t{p['start']}\t{p['end']}\t{p['summit']}\t{p['center']}\t"
                     f"{p['chrom']}:{p['ws'] + 1}-{p['ws'] + input_len}\t{len(p['variants'])}\t{p['n_indel']}\t"
                     f"{p['n_multi']}\t{p['n_refmis']}\t{len(p['keys'])}\t{p['n_ambiguous']}\t{p['n_missing']}\n")

    # one row per unique haplotype of each peak (hap_id 0 = hg38 REF)
    with gzip.open(os.path.join(chunk_dir, "haplotypes.tsv.gz"), "wt") as fh:
        fh.write("\t".join(["peak_id", "hap_id", "n_donor_haplotypes", "n_alt", "alt_alleles"] + cols) + "\n")
        for p in peaks:
            for k, key in enumerate(p["keys"]):
                g = p["offset"] + k
                alts = ",".join(f"{p['variants'][v][0]}{p['variants'][v][1]}>{p['variants'][v][2]}"
                                for v in key) or "REF"
                fh.write("\t".join([p["name"], str(k), str(p["n_carriers"][k]), str(len(key)), alts]
                                   + [f"{counts[c][g]:.4f}" for c in cols]) + "\n")

    # each donor's two hap_ids per peak ("hap1|hap2"), to link donors to haplotypes.tsv.gz
    with gzip.open(os.path.join(chunk_dir, "donor_hap_ids.tsv.gz"), "wt") as fh:
        fh.write("peak_id\t" + "\t".join(samples) + "\n")
        for p in peaks:
            fh.write(p["name"] + "\t" + "\t".join(f"{a}|{b}" for a, b in p["hap_ids"]) + "\n")

    # peaks x donors matrices: diploid predicted counts = mean of the two haplotypes
    for c in matrix_cols:
        with gzip.open(os.path.join(chunk_dir, f"diploid_{c}.tsv.gz"), "wt") as fh:
            fh.write("peak_id\t" + "\t".join(samples) + "\n")
            for p in peaks:
                g = p["offset"] + p["hap_ids"]
                dip = (counts[c][g[:, 0]] + counts[c][g[:, 1]]) / 2
                fh.write(p["name"] + "\t" + "\t".join(f"{x:.2f}" for x in dip) + "\n")
    print(f"Finalized {len(peaks)} peaks x {len(samples)} donors, {len(matrix_cols)} matrices in {time.time() - t:.1f} s")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "predict":
        predict(sys.argv[2:])
    elif mode == "finalize":
        finalize(sys.argv[2:])
    else:
        sys.exit(f"unknown mode {mode}")


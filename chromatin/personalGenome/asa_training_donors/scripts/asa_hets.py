"""
Step 1b, part 1: heterozygous SNVs and predicted haplotype counts per training donor, from the Step 1 outputs.

The het sites are taken from what the model saw: a donor's two haplotypes (donor_hap_ids.tsv.gz, hap1 = first
allele of the phased GT) and the ALT alleles each haplotype carries (haplotypes.tsv.gz). A het = an ALT allele
carried by exactly one of the two haplotypes. Only hets inside the model's central 1,000 bp output region
(summit - 499 .. summit + 500, 1-based) are kept, because that is the region the predicted counts describe.

Usage (plain Python 3, no numpy needed):
  python asa_hets.py <step1_merged_dir> <cell_type> <id_map.tsv> <out_dir>
    id_map.tsv : 2 columns, BAM donor name -> VCF sample ID (only these donors are processed)

Writes to <out_dir>:
  <donor>.sites.tsv      chrom, pos, ref, alt, hap1_alt (1 = hap1 carries ALT), peak_id, core_start, core_end
  <donor>.positions.txt  chrom, pos (unique, for samtools mpileup -l)
  pred.tsv               donor, peak_id, n_het_core, hap1_id, hap2_id, <cond>_hap1, <cond>_hap2 per condition
                         (fold means), then <cond>_fold<N>_hap1, <cond>_fold<N>_hap2 per fold
  DONE
"""
import gzip
import os
import sys

CORE_HALF = 500        # model output = central 1,000 bp


def main(merged, ct, id_map, out):
    os.makedirs(out, exist_ok=True)
    donors = [l.split()[:2] for l in open(id_map) if len(l.split()) >= 2 and not l.startswith("#")]

    # --- peak cores (center column = 0-based summit position) ---
    core = {}
    with open(os.path.join(merged, "peak_info.tsv")) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ic, ich = h.index("center"), h.index("chrom")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            c = int(f[ic])
            core[f[0]] = (f[ich], c - CORE_HALF + 1, c + CORE_HALF)

    # --- each donor's hap1|hap2 per peak; keep peaks where the two haplotypes differ ---
    with gzip.open(os.path.join(merged, "donor_hap_ids.tsv.gz"), "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        missing = [v for _, v in donors if v not in h]
        if missing:
            sys.exit(f"ERROR VCF sample IDs not in donor_hap_ids.tsv.gz header: {missing}")
        cols = [(b, h.index(v)) for b, v in donors]
        het_peaks = {b: [] for b, _ in donors}          # donor -> [(peak, hap1, hap2)]
        needed = set()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            for b, i in cols:
                h1, h2 = f[i].split("|")
                if h1 != h2:
                    het_peaks[b].append((f[0], h1, h2))
                    needed.add((f[0], h1))
                    needed.add((f[0], h2))

    # --- ALT alleles and predicted counts (fold means, then each fold) of the haplotypes we need ---
    hap = {}
    with gzip.open(os.path.join(merged, "haplotypes.tsv.gz"), "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        mean_cols = [(c[:-5], h.index(c)) for c in h if c.endswith("_mean")]    # e.g. PBS_mean -> PBS
        mean_cols += [(c, h.index(c)) for c in h if "_fold" in c]                 # e.g. PBS_fold0
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if (f[0], f[1]) in needed:
                alts = set() if f[4] == "REF" else set(f[4].split(","))
                hap[(f[0], f[1])] = (alts, [f[i] for _, i in mean_cols])
    conds = [c for c, _ in mean_cols]

    # --- per donor: het sites in the core + prediction table ---
    with open(os.path.join(out, "pred.tsv"), "w") as fp:
        fp.write("donor\tpeak_id\tn_het_core\thap1_id\thap2_id\t"
                 + "\t".join(f"{c}_hap1\t{c}_hap2" for c in conds) + "\n")
        for b, _ in donors:
            sites = []
            for peak, h1, h2 in het_peaks[b]:
                chrom, cs, ce = core[peak]
                a1, p1 = hap[(peak, h1)]
                a2, p2 = hap[(peak, h2)]
                n = 0
                for allele in a1 ^ a2:                    # "25020562C>T"
                    pos, ra = int(allele[:-3]), allele[-3:]
                    if cs <= pos <= ce:
                        sites.append((chrom, pos, ra[0], ra[2], int(allele in a1), peak, cs, ce))
                        n += 1
                if n:
                    fp.write(f"{b}\t{peak}\t{n}\t{h1}\t{h2}\t"
                             + "\t".join(f"{x}\t{y}" for x, y in zip(p1, p2)) + "\n")
            sites.sort(key=lambda s: (s[0], s[1]))
            with open(os.path.join(out, f"{b}.sites.tsv"), "w") as fs:
                fs.write("chrom\tpos\tref\talt\thap1_alt\tpeak_id\tcore_start\tcore_end\n")
                for s in sites:
                    fs.write("\t".join(map(str, s)) + "\n")
            with open(os.path.join(out, f"{b}.positions.txt"), "w") as fs:
                for chrom, pos in sorted({(s[0], s[1]) for s in sites}):
                    fs.write(f"{chrom}\t{pos}\n")
            print(f"{ct} {b}: {len(het_peaks[b])} peaks with two different haplotypes; "
                  f"{len(sites)} het sites in the central 1 kb")
    open(os.path.join(out, "DONE"), "w").close()


if __name__ == "__main__":
    main(*sys.argv[1:5])

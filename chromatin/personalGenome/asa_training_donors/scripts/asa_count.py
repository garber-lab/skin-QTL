"""
Step 1b, part 2: allele and haplotype read counts for one BAM, from `samtools mpileup --output-QNAME` on stdin.

  samtools mpileup -l <donor>.positions.txt -q 30 -Q 20 -B -d 100000 --output-QNAME <bam> \
      | python asa_count.py <donor>.sites.tsv <out_prefix> <context> <donor>

mpileup already drops duplicates, secondary/QC-fail reads, MAPQ < 30 and bases with quality < 20, and counts
the overlapping part of a read pair once. Here each read (pair) is assigned to a haplotype through the phased
het it covers; a pair covering several hets of the same peak is counted once per peak (pairs whose alleles
point to different haplotypes are counted as 'conflict' and left out).

Writes:
  <out_prefix>.sites.tsv   per het site: read counts for REF, ALT and other bases
  <out_prefix>.peaks.tsv   per peak: read pairs on hap1 / hap2 / conflicting
"""
import sys


def parse_bases(bases):
    """mpileup base column -> one entry per read: 'R' (= reference), a base letter, or '' (deletion etc.)."""
    out, i, n = [], 0, len(bases)
    while i < n:
        ch = bases[i]
        if ch == "^":                          # read start, next char = MAPQ
            i += 2
            continue
        if ch == "$":                          # read end
            i += 1
            continue
        if ch in "+-":                         # indel after the previous base: +3ACG
            j = i + 1
            while j < n and bases[j].isdigit():
                j += 1
            i = j + int(bases[i + 1:j])
            continue
        if ch in ".,":
            out.append("R")
        elif ch.upper() in "ACGTN":
            out.append(ch.upper())
        else:                                  # * # > < : deletion / reference skip
            out.append("")
        i += 1
    return out


def main(sites_tsv, prefix, context, donor):
    # sites: (chrom, pos) -> REF, ALT, hap1_alt, [peaks]; peaks -> (chrom, core_end)
    sites, peak_end = {}, {}
    with open(sites_tsv) as fh:
        fh.readline()
        for line in fh:
            chrom, pos, ref, alt, h1alt, peak, cs, ce = line.rstrip("\n").split("\t")
            key = (chrom, int(pos))
            if key not in sites:
                sites[key] = [ref, alt, int(h1alt), []]
            sites[key][3].append(peak)
            peak_end[peak] = (chrom, int(ce))

    site_counts = {}                           # (chrom, pos) -> [ref, alt, other]
    peak_counts = {}                           # peak -> [n_sites_covered, hap1, hap2, conflict]
    active = {}                                # peak -> {qname: set(haps)}, peaks not yet passed

    def flush(peak):
        reads = active.pop(peak)
        c = peak_counts.setdefault(peak, [0, 0, 0, 0])
        for haps in reads.values():
            if len(haps) == 1:
                c[1 if 1 in haps else 2] += 1
            else:
                c[3] += 1

    cur_chrom = None
    for line in sys.stdin:
        f = line.rstrip("\n").split("\t")
        chrom, pos = f[0], int(f[1])
        if chrom != cur_chrom:
            for p in list(active):
                flush(p)
            cur_chrom = chrom
        else:
            for p in [p for p in active if peak_end[p][1] < pos]:
                flush(p)
        s = sites.get((chrom, pos))
        # depth 0 = all reads at the site were filtered (MAPQ/BQ); samtools then prints '*' in the base,
        # quality and name columns
        if s is None or len(f) < 7 or f[3] == "0":
            continue
        ref, alt, h1alt, peaks = s
        calls = parse_bases(f[4])
        qnames = f[6].split(",") if f[6] not in ("", "*") else []
        if len(qnames) != len(calls):
            sys.exit(f"ERROR {chrom}:{pos}: {len(calls)} bases but {len(qnames)} read names")
        sc = site_counts.setdefault((chrom, pos), [0, 0, 0])
        for call, q in zip(calls, qnames):
            if call == "R" or call == ref:
                sc[0] += 1
                h = 2 if h1alt else 1              # read carries REF -> the haplotype with REF
            elif call == alt:
                sc[1] += 1
                h = 1 if h1alt else 2
            else:
                if call:
                    sc[2] += 1
                continue
            for p in peaks:
                active.setdefault(p, {}).setdefault(q, set()).add(h)
        for p in peaks:
            pc = peak_counts.setdefault(p, [0, 0, 0, 0])
            pc[0] += 1
    for p in list(active):
        flush(p)

    with open(f"{prefix}.sites.tsv", "w") as fh:
        fh.write("context\tdonor\tchrom\tpos\tref\talt\thap1_alt\tpeak_ids\tref_reads\talt_reads\tother_reads\n")
        for key in sorted(site_counts):
            ref, alt, h1alt, peaks = sites[key]
            r, a, o = site_counts[key]
            fh.write(f"{context}\t{donor}\t{key[0]}\t{key[1]}\t{ref}\t{alt}\t{h1alt}\t{','.join(peaks)}\t{r}\t{a}\t{o}\n")
    with open(f"{prefix}.peaks.tsv", "w") as fh:
        fh.write("context\tdonor\tpeak_id\tn_sites_covered\thap1_reads\thap2_reads\tconflict_reads\n")
        for p in sorted(peak_counts):
            fh.write(f"{context}\t{donor}\t{p}\t" + "\t".join(map(str, peak_counts[p])) + "\n")
    n_r = sum(c[0] for c in site_counts.values())
    n_a = sum(c[1] for c in site_counts.values())
    print(f"{context} {donor}: {len(site_counts)} covered het sites of {len(sites)}; "
          f"{len(peak_counts)} peaks; REF fraction {n_r / max(n_r + n_a, 1):.3f}")


if __name__ == "__main__":
    main(*sys.argv[1:5])

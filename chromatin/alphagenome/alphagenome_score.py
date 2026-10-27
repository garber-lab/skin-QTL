#!/usr/bin/env python
"""Score a list of variants with the hosted AlphaGenome API.

Input: tab-separated file with columns
    variant_id  CHROM  POS  REF  ALT
CHROM must start with 'chr' (e.g. chr1). POS is 1-based, hg38.

Output: tab-separated table (gzipped if the name ends in .gz), one row per
variant x scorer x track (x gene, for gene-level scorers).

Usage:
    python alphagenome_score.py --input variants.tsv --output scores.tsv.gz

Inside an LSF job array, the script automatically processes only its share of
the variants (every Nth row), using $LSB_JOBINDEX and --n-chunks.
"""

import argparse
import os
import sys
import time

import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', required=True, help='Variant table (TSV).')
    p.add_argument('--output', required=True, help='Output table (.tsv or .tsv.gz).')
    p.add_argument('--sequence-length', default='1MB',
                   choices=['16KB', '100KB', '500KB', '1MB'],
                   help='DNA window centred on each variant (default: 1MB).')
    p.add_argument('--scorers', default='rna_seq,atac,dnase,chip_histone',
                   help='Comma-separated scorer names (case-insensitive), or "all".')
    p.add_argument('--n-chunks', type=int, default=1,
                   help='Total number of job-array tasks splitting the input.')
    p.add_argument('--chunk-index', type=int,
                   default=int(os.environ.get('LSB_JOBINDEX', '0') or 0),
                   help='1-based task index (defaults to $LSB_JOBINDEX).')
    p.add_argument('--max-retries', type=int, default=3,
                   help='Retries per variant on API/network errors.')
    return p.parse_args()


def main():
    args = parse_args()

    api_key = os.environ.get('ALPHAGENOME_API_KEY')
    if not api_key:
        sys.exit('ERROR: set the ALPHAGENOME_API_KEY environment variable.')

    # ---- Read and check the variant table ----
    vcf = pd.read_csv(args.input, sep='\t')
    required = ['variant_id', 'CHROM', 'POS', 'REF', 'ALT']
    missing = [c for c in required if c not in vcf.columns]
    if missing:
        sys.exit(f'ERROR: input is missing columns: {missing}')

    # ---- Keep only this task's share of rows when run as a job array ----
    if args.n_chunks > 1:
        if not 1 <= args.chunk_index <= args.n_chunks:
            sys.exit(f'ERROR: chunk index {args.chunk_index} not in 1..{args.n_chunks}')
        vcf = vcf.iloc[args.chunk_index - 1::args.n_chunks]
    print(f'Scoring {len(vcf)} variants', flush=True)

    # ---- Choose scorers ----
    all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    if args.scorers.strip().lower() == 'all':
        selected = list(all_scorers.values())
    else:
        wanted = {s.strip().lower() for s in args.scorers.split(',')}
        selected = [all_scorers[k] for k in all_scorers if k.lower() in wanted]
        unknown = wanted - {k.lower() for k in all_scorers}
        if unknown or not selected:
            sys.exit(f'ERROR: unknown scorers {sorted(unknown)}. '
                     f'Available: {[k.lower() for k in all_scorers]}')
    print('Scorers:', [str(s) for s in selected], flush=True)

    seq_len = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
        f'SEQUENCE_LENGTH_{args.sequence_length}']
    organism = dna_client.Organism.HOMO_SAPIENS

    model = dna_client.create(api_key)

    # ---- Score each variant, retrying transient failures ----
    results, failed = [], []
    for i, row in enumerate(vcf.itertuples(index=False), start=1):
        variant = genome.Variant(
            chromosome=str(row.CHROM),
            position=int(row.POS),
            reference_bases=row.REF,
            alternate_bases=row.ALT,
            name=str(row.variant_id),
        )
        interval = variant.reference_interval.resize(seq_len)

        for attempt in range(1, args.max_retries + 1):
            try:
                scores = model.score_variant(
                    interval=interval,
                    variant=variant,
                    variant_scorers=selected,
                    organism=organism,
                )
                results.append(scores)
                break
            except Exception as e:  # network hiccups, rate limits, bad input
                print(f'  {row.variant_id}: attempt {attempt} failed: {e}', flush=True)
                if attempt == args.max_retries:
                    failed.append((row.variant_id, str(e)))
                else:
                    time.sleep(10 * attempt)

        if i % 10 == 0 or i == len(vcf):
            print(f'{i}/{len(vcf)} done ({len(failed)} failed)', flush=True)

    # ---- Write outputs ----
    if results:
        df = variant_scorers.tidy_scores(results)
        df.to_csv(args.output, sep='\t', index=False)  # .gz => gzipped
        print(f'Wrote {len(df)} rows to {args.output}', flush=True)
    else:
        print('No variants scored successfully; no output table written.', flush=True)

    if failed:
        fail_path = args.output + '.failed.tsv'
        pd.DataFrame(failed, columns=['variant_id', 'error']).to_csv(
            fail_path, sep='\t', index=False)
        print(f'{len(failed)} variants failed; listed in {fail_path}', flush=True)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Construct site-by-SV ALT allele frequency matrices from a biallelic SV VCF.

Install: python -m pip install pandas openpyxl numpy
Run: python 01_make_site_AF_v3.py --vcf sv_envAssoc.filtered.noHighMiss4.recalc.vcf.gz --samples sample_site.xlsx --outdir AF_output

Site AF is computed as ALT allele copies divided by total *called* allele copies.
Unique keys comprise CHROM, POS and VCF ID only. Colliding keys are rejected
rather than silently merging variants. REF and ALT remain in metadata.
"""
import argparse
import gzip
import re
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd


def open_vcf(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'rt')


def parse_info(info, tag):
    m = re.search(r'(?:^|;)' + re.escape(tag) + r'=([^;]+)', info)
    return m.group(1) if m else ''


def parse_gt(sample_text, gt_index):
    chunks = sample_text.split(':')
    if gt_index is None or gt_index >= len(chunks):
        return None
    gt = chunks[gt_index]
    alleles = re.split(r'[/|]', gt)
    if not alleles or any(a not in ('0', '1') for a in alleles):
        return None
    return sum(a == '1' for a in alleles), len(alleles)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', default='sv_envAssoc.filtered.noHighMiss4.recalc.vcf.gz')
    parser.add_argument('--samples', default='sample_site.xlsx')
    parser.add_argument('--outdir', default='AF_output')
    parser.add_argument('--min-called-site', type=int, default=1,
                        help='Minimum number of called individuals per site for AF; set to 2/3 for sensitivity checks.')
    args = parser.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    meta = pd.read_excel(args.samples, dtype={'sample_id': str, 'site_id': str})
    required = {'sample_id', 'site_id'}
    if not required.issubset(meta.columns):
        raise ValueError(f'Missing metadata columns: {sorted(required - set(meta.columns))}')
    meta = meta.dropna(subset=['sample_id', 'site_id']).copy()
    for c in ['sample_id', 'site_id']:
        meta[c] = meta[c].str.strip()
    if meta['sample_id'].duplicated().any():
        raise ValueError('Duplicate sample_id in sample_site.xlsx: ' + ', '.join(meta.loc[meta['sample_id'].duplicated(), 'sample_id'].tolist()))
    if not len(meta):
        raise ValueError('Empty sample metadata')
    samples = meta['sample_id'].tolist()
    sites = sorted(meta['site_id'].unique())
    site_to_indices = {s: np.flatnonzero(meta['site_id'].to_numpy() == s) for s in sites}
    meta.groupby('site_id', sort=True).size().rename('n_samples').reset_index().to_csv(out / '00_site_sample_counts.tsv', sep='\t', index=False)

    af_data = {s: [] for s in sites}
    called_data = {s: [] for s in sites}
    an_data = {s: [] for s in sites}
    records = []
    used_keys = set()
    duplicate_raw_id = defaultdict(int)
    skipped_multiallelic = 0
    skipped_no_gt = 0
    ploidies = set()
    n_variants = 0
    ordered_idx = None
    with open_vcf(args.vcf) as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header = line.rstrip('\n').split('\t')
                vcf_samples = header[9:]
                if len(set(vcf_samples)) != len(vcf_samples):
                    raise ValueError('Duplicate sample IDs in VCF header')
                sidx = {s: i for i, s in enumerate(vcf_samples)}
                missing = [s for s in samples if s not in sidx]
                if missing:
                    raise ValueError('Metadata sample_id absent from VCF: ' + ', '.join(missing))
                ordered_idx = [sidx[s] for s in samples]
                continue
            if not line.strip() or line.startswith('#'):
                continue
            if ordered_idx is None:
                raise ValueError('No #CHROM VCF header found')
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 10:
                raise ValueError(f'VCF record with fewer than 10 columns: {line[:120]}')
            chrom, pos_text, raw_id, ref, alt = fields[:5]
            if ',' in alt:
                skipped_multiallelic += 1
                continue
            fmt = fields[8].split(':')
            if 'GT' not in fmt:
                skipped_no_gt += 1
                continue
            gtidx = fmt.index('GT')
            pos = int(pos_text)
            # Short identifier requested by the user. Check uniqueness strictly.
            key = f'{chrom}|{pos}|{raw_id}'
            if key in used_keys:
                raise ValueError(
                    f'Duplicate CHROM+POS+ID key: {key}. '
                    'These VCF rows cannot be uniquely identified by the requested short key. '
                    'Inspect those records before continuing.'
                )
            used_keys.add(key)
            duplicate_raw_id[raw_id] += 1
            parsed = [parse_gt(fields[9 + i], gtidx) for i in ordered_idx]
            for g in parsed:
                if g is not None:
                    ploidies.add(g[1])
            end_text = parse_info(fields[7], 'END')
            svtype = parse_info(fields[7], 'SVTYPE')
            end = int(end_text) if end_text.isdigit() else (pos + max(len(ref) - 1, 0) if not svtype.startswith('INS') else pos)
            rec = {'variant_key': key, 'analysis_sv_id': raw_id, 'chromosome': chrom,
                   'POS': pos, 'END': end, 'SVTYPE': svtype, 'REF': ref, 'ALT': alt}
            records.append(rec)
            for site, idx in site_to_indices.items():
                valid = [parsed[i] for i in idx if parsed[i] is not None]
                called_n = len(valid)
                called_alleles = sum(v[1] for v in valid)
                alt_copies = sum(v[0] for v in valid)
                af = alt_copies / called_alleles if called_n >= args.min_called_site and called_alleles > 0 else np.nan
                af_data[site].append(af)
                called_data[site].append(called_n)
                an_data[site].append(called_alleles)
            n_variants += 1
            if n_variants % 10000 == 0:
                print(f'Processed {n_variants:,} variants...', flush=True)

    keys = [r['variant_key'] for r in records]
    if not keys:
        raise ValueError('No valid biallelic variants with GT found')
    pd.DataFrame(records).to_csv(out / '04_SV_metadata.tsv', sep='\t', index=False)
    for filename, source in [('02_site_SV_AF.tsv', af_data),
                             ('03_site_SV_calledN.tsv', called_data),
                             ('03b_site_SV_calledAN.tsv', an_data)]:
        mat = pd.DataFrame.from_dict(source, orient='index', columns=keys)
        mat.index.name = 'site_id'
        mat.to_csv(out / filename, sep='\t', na_rep='NA', float_format='%.9g')
    dup_ids = pd.DataFrame([{'analysis_sv_id': rid, 'vcf_record_count': count}
                            for rid, count in duplicate_raw_id.items() if count > 1])
    pd.DataFrame(dup_ids, columns=['analysis_sv_id', 'vcf_record_count']).to_csv(out / '05_duplicate_raw_IDs.tsv', sep='\t', index=False)
    print(f'Done: {n_variants:,} retained variants, {len(samples)} individuals, {len(sites)} sites')
    print(f'Duplicate raw IDs: {len(dup_ids)}; skipped multiallelic: {skipped_multiallelic}; skipped records without GT: {skipped_no_gt}')
    print(f'Called ploidies: {sorted(ploidies)}. X chromosome AF denominators account for observed GT ploidy.')
    print(f'Output directory: {out.resolve()}')

if __name__ == '__main__':
    main()

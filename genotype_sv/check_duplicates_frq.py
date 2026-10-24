#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check duplicates in vcftools --freq output (frq). 
Usage:
  python check_duplicates_frq.py china.frq
Outputs:
  - prints summary counts
  - writes a deduplicated file <input>.dedup.frq when --dedup is provided
Dedup strategy: keep rows by (CHROM,POS,ALLELE) uniqueness; if multiple identical (CHROM,POS,ALLELE),
we take the first occurrence (you can modify to average).
"""
import sys
import re
from collections import defaultdict

def parse_line_into_alleles(line):
    # split into 5 parts max
    parts = line.rstrip("\n").split("\t", 4)
    if len(parts) < 5:
        return None
    chrom, pos, n_alleles, n_chr, allele_field = parts
    pos = int(pos)
    # allele_field may contain multiple allele:FREQ separated by whitespace or tabs
    # split by whitespace
    allele_tokens = re.split(r"\s+", allele_field.strip())
    alleles = []
    for tok in allele_tokens:
        if ':' in tok:
            allele, freq = tok.split(":", 1)
            try:
                af = float(freq)
            except:
                af = None
            alleles.append((allele, af, tok))
    return chrom, pos, alleles, line

def check_file(fp, dedup=False):
    seen_pos = defaultdict(int)
    seen_triplet = {}
    duplicates_examples = []
    total_lines = 0
    total_alleles = 0
    out_lines = []
    with open(fp, 'r', encoding='utf-8') as f:
        header = f.readline()
        for line in f:
            total_lines += 1
            parsed = parse_line_into_alleles(line)
            if parsed is None:
                continue
            chrom, pos, alleles, raw = parsed
            seen_pos[(chrom,pos)] += 1
            for allele, af, tok in alleles:
                total_alleles += 1
                key = (chrom,pos,allele)
                if key in seen_triplet:
                    # duplicate allele occurrence
                    duplicates_examples.append((key, seen_triplet[key], af))
                else:
                    seen_triplet[key] = af
                    if dedup:
                        # reconstruct a minimal line: we keep original raw for first occurrence
                        out_lines.append((chrom,pos,allele,af,raw))
    print(f"File: {fp}")
    print(f"  Total input lines (excluding header): {total_lines}")
    print(f"  Total allele tokens parsed: {total_alleles}")
    dup_pos_count = sum(1 for k,v in seen_pos.items() if v>1)
    print(f"  Number of (CHROM,POS) with >1 lines: {dup_pos_count}")
    print(f"  Number of duplicate (CHROM,POS,ALLELE) occurrences: {len(duplicates_examples)}")
    if len(duplicates_examples)>0:
        print("  Example duplicate entries (CHROM,POS,ALLELE), existing_af, new_af:")
        for ex in duplicates_examples[:10]:
            print("   ", ex)
    if dedup:
        outname = fp + ".dedup.frq"
        # write header + unique lines (we stored raw lines for first occurrence)
        with open(outname, 'w', encoding='utf-8') as out:
            out.write(header)
            # out_lines contains tuples with raw lines; they are in first-seen order
            for _,_,_,_,raw in out_lines:
                out.write(raw)
        print("  Wrote deduplicated file:", outname)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python check_duplicates_frq.py china.frq [--dedup]")
        sys.exit(1)
    fp = sys.argv[1]
    dedup_flag = "--dedup" in sys.argv
    check_file(fp, dedup=dedup_flag)


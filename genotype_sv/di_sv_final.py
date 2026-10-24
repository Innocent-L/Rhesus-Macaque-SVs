#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import re

def parse_vcftools_frq(file):
    data = {}
    with open(file, encoding="utf-8") as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split()
            chrom, pos = parts[0], parts[1]
            allele_fields = parts[4:]
            alt_af = None
            for field in allele_fields:
                allele, af = field.split(':')
                # input
                if not re.match(r'^[ACGTN]$', allele):
                    alt_af = float(af)
                    break
            if alt_af is not None:
                data[(chrom, int(pos))] = alt_af
    return data


china_file = "china.frq"
india_file = "india.frq"

print("Parsing China freq...")
china_af = parse_vcftools_frq(china_file)
print("Parsing India freq...")
india_af = parse_vcftools_frq(india_file)

common_keys = set(china_af.keys()) & set(india_af.keys())
print(f"share {len(common_keys)} DI_SV ")

# --- sort ---
def chrom_order(chrom):
    m = re.search(r'(\d+)$', chrom)
    return int(m.group(1)) if m else 1e9  

sorted_keys = sorted(common_keys, key=lambda x: (chrom_order(x[0]), x[1]))

with open("DI_SV_siteLevel.sorted.tsv", "w", encoding="utf-8") as outf:
    outf.write("CHROM\tPOS\tAF_china\tAF_india\tDI_SV\n")
    for k in sorted_keys:
        af_china = china_af[k]
        af_india = india_af[k]
        di = abs(af_china - af_india)
        outf.write(f"{k[0]}\t{k[1]}\t{af_china}\t{af_india}\t{di}\n")

print("Done. Output -> DI_SV_siteLevel.sorted.tsv")


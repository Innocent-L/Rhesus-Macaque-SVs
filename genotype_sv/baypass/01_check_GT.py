#!/usr/bin/env python3

import re

VCF = "all_samples_rename.AF_filter.vcf"

REMOVE = {
    # long-read individuals
    "cx_2",
    "fj",
    "hb_1",
    "xz_1",
    "zm_2",
    "yd_10",

    # Indian population individuals
    "yd_117",
    "yd_141",
    "yd_144",
    "yd_105",
    "yd_121",
    "yd_275",
    "yd_153",
    "yd_147",
    "yd_125",
    "yd_133",
}

VALID_GT = {
    "0/0", "0/1", "1/0", "1/1",
    "0|0", "0|1", "1|0", "1|1",
    ".", "./.", ".|."
}

samples = None
keep_indices = None

bad = []

with open(VCF, "r") as f:

    for line in f:

        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):

            x = line.rstrip("\n").split("\t")

            samples = x[9:]

            keep_indices = [
                i for i, s in enumerate(samples)
                if s not in REMOVE
            ]

            kept_samples = [
                samples[i]
                for i in keep_indices
            ]

            print("Retained samples:", len(kept_samples))

            continue

        if line.startswith("#"):
            continue

        x = line.rstrip("\n").split("\t")

        chrom = x[0]
        pos = x[1]
        vid = x[2]

        fmt = x[8].split(":")

        if "GT" not in fmt:

            bad.append(
                (chrom, pos, vid, "NO_GT_FIELD", "")
            )

            continue

        gt_index = fmt.index("GT")

        sample_fields = x[9:]

        for i in keep_indices:

            sample = samples[i]

            values = sample_fields[i].split(":")

            if gt_index >= len(values):
                gt = ""
            else:
                gt = values[gt_index]

            if gt not in VALID_GT:

                bad.append(
                    (
                        chrom,
                        pos,
                        vid,
                        sample,
                        repr(gt)
                    )
                )


with open("invalid_GT_china46.tsv", "w") as out:

    out.write(
        "CHROM\tPOS\tID\tSAMPLE\tRAW_GT\n"
    )

    for row in bad:

        out.write(
            "\t".join(map(str, row))
            + "\n"
        )


print("Invalid GT records:", len(bad))

if bad:

    print("\nFirst 20 invalid GT records:")

    for row in bad[:20]:
        print(*row, sep="\t")

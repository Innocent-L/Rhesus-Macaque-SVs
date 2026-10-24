#!/usr/bin/env python3

import os
import re


############################################################
# Input / output
############################################################

VCF = "all_samples_rename.AF_filter.vcf"

OUTDIR = "BayPass_China46"

os.makedirs(OUTDIR, exist_ok=True)


GENO_OUT = os.path.join(
    OUTDIR,
    "china46_allpoly.geno"
)

MARKER_OUT = os.path.join(
    OUTDIR,
    "china46_allpoly.marker.tsv"
)

QC_OUT = os.path.join(
    OUTDIR,
    "china46_allpoly.qc.tsv"
)

INVALID_GT_OUT = os.path.join(
    OUTDIR,
    "china46_invalid_GT.tsv"
)


############################################################
# Individuals to remove
############################################################

REMOVE = {

    # six long-read / assembly individuals
    "cx_2",
    "fj",
    "hb_1",
    "xz_1",
    "zm_2",
    "yd_10",

    # ten Indian population individuals
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


############################################################
# Population order
#
# IMPORTANT:
# this order determines BayPass columns
############################################################

POP_ORDER = [
    "cx",
    "fj",
    "hb",
    "hn",
    "zm"
]


EXPECTED_N = {
    "cx": 12,
    "fj": 10,
    "hb": 7,
    "hn": 8,
    "zm": 9
}


############################################################
# Valid diploid biallelic GT
############################################################

VALID_MISSING = {
    ".",
    "./.",
    ".|.",
    ""
}


VALID_CALLED = {
    "0/0",
    "0/1",
    "1/0",
    "1/1",
    "0|0",
    "0|1",
    "1|0",
    "1|1"
}


############################################################
# QC counters
############################################################

qc = {
    "VCF_records_total": 0,
    "multiallelic_skipped": 0,
    "no_GT_field_skipped": 0,
    "all_missing_skipped": 0,
    "monomorphic_in_China_skipped": 0,
    "retained_polymorphic_SVs": 0,
    "invalid_GT_calls_set_to_missing": 0
}


############################################################
# Storage
############################################################

samples = None
keep_indices = None
sample_pop = {}

invalid_gt_records = []


############################################################
# Outputs
############################################################

geno_out = open(
    GENO_OUT,
    "w"
)

marker_out = open(
    MARKER_OUT,
    "w"
)


header = [
    "MRK",
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "FILTER",
    "END",
    "SVTYPE",
    "SVLEN",
    "TOTAL_REF",
    "TOTAL_ALT",
    "TOTAL_AN",
    "GLOBAL_ALT_AF"
]


for p in POP_ORDER:

    header.extend([
        f"{p}_REF",
        f"{p}_ALT",
        f"{p}_AN",
        f"{p}_AF",
        f"{p}_N_called"
    ])


marker_out.write(
    "\t".join(header)
    + "\n"
)


############################################################
# Parse VCF directly
############################################################

with open(VCF, "r") as f:

    mrk = 0

    for line in f:

        ####################################################
        # Header
        ####################################################

        if line.startswith("##"):
            continue


        if line.startswith("#CHROM"):

            x = line.rstrip("\n").split("\t")

            samples = x[9:]


            ################################################
            # Keep China short-read individuals only
            ################################################

            keep_indices = [
                i
                for i, sample in enumerate(samples)
                if sample not in REMOVE
            ]


            kept_samples = [
                samples[i]
                for i in keep_indices
            ]


            if len(kept_samples) != 46:

                raise RuntimeError(
                    f"Expected 46 retained samples, "
                    f"found {len(kept_samples)}"
                )


            ################################################
            # Population assignment
            ################################################

            pop_n = {
                p: 0
                for p in POP_ORDER
            }


            for sample in kept_samples:

                pop = sample.split("_")[0]

                if pop not in POP_ORDER:

                    raise RuntimeError(
                        f"Unexpected retained population: "
                        f"{sample} -> {pop}"
                    )

                sample_pop[sample] = pop
                pop_n[pop] += 1


            if pop_n != EXPECTED_N:

                raise RuntimeError(
                    f"Population counts incorrect.\n"
                    f"Observed: {pop_n}\n"
                    f"Expected: {EXPECTED_N}"
                )


            print(
                "Retained population sizes:",
                pop_n
            )

            continue


        if line.startswith("#"):
            continue


        ####################################################
        # Variant record
        ####################################################

        qc["VCF_records_total"] += 1

        x = line.rstrip("\n").split("\t")


        if len(x) < 10:
            continue


        chrom = x[0]
        pos = x[1]
        vid = x[2]

        ref = x[3]
        alt = x[4]

        filt = x[6]
        info_string = x[7]


        ####################################################
        # BayPass is biallelic
        ####################################################

        if "," in alt:

            qc["multiallelic_skipped"] += 1
            continue


        ####################################################
        # INFO only for metadata
        ####################################################

        info = {}

        for z in info_string.split(";"):

            if "=" in z:

                key, value = z.split(
                    "=",
                    1
                )

                info[key] = value


        end = info.get(
            "END",
            "."
        )

        svtype = info.get(
            "SVTYPE",
            "."
        )

        svlen = info.get(
            "SVLEN",
            "."
        )


        ####################################################
        # Locate GT
        ####################################################

        fmt = x[8].split(":")


        if "GT" not in fmt:

            qc["no_GT_field_skipped"] += 1
            continue


        gt_index = fmt.index("GT")

        sample_fields = x[9:]


        ####################################################
        # Per-population allele counts
        ####################################################

        ref_count = {
            p: 0
            for p in POP_ORDER
        }

        alt_count = {
            p: 0
            for p in POP_ORDER
        }

        n_called = {
            p: 0
            for p in POP_ORDER
        }


        ####################################################
        # Iterate retained samples
        ####################################################

        for i in keep_indices:

            sample = samples[i]

            pop = sample_pop[sample]


            if i >= len(sample_fields):

                gt = ""

            else:

                values = sample_fields[i].split(":")

                if gt_index >= len(values):

                    gt = ""

                else:

                    gt = values[gt_index]


            ################################################
            # Missing
            ################################################

            if gt in VALID_MISSING:
                continue


            ################################################
            # Invalid syntactic GT
            #
            # Do NOT infer genotype.
            # Treat it as missing and record it.
            ################################################

            if gt not in VALID_CALLED:

                qc[
                    "invalid_GT_calls_set_to_missing"
                ] += 1


                invalid_gt_records.append([
                    chrom,
                    pos,
                    vid,
                    sample,
                    gt
                ])

                continue


            ################################################
            # Valid genotype
            ################################################

            alleles = re.split(
                r"[/|]",
                gt
            )


            ref_count[pop] += (
                alleles.count("0")
            )

            alt_count[pop] += (
                alleles.count("1")
            )

            n_called[pop] += 1


        ####################################################
        # Total allele counts
        ####################################################

        total_ref = sum(
            ref_count.values()
        )

        total_alt = sum(
            alt_count.values()
        )

        total_an = (
            total_ref
            + total_alt
        )


        ####################################################
        # No usable genotype
        ####################################################

        if total_an == 0:

            qc["all_missing_skipped"] += 1
            continue


        ####################################################
        # Remove variants fixed across all Chinese samples
        #
        # REF fixed:
        # ALT = 0
        #
        # ALT fixed:
        # REF = 0
        #
        # Neither contributes to differentiation.
        ####################################################

        if (
            total_ref == 0
            or
            total_alt == 0
        ):

            qc[
                "monomorphic_in_China_skipped"
            ] += 1

            continue


        ####################################################
        # Global ALT AF
        #
        # Descriptive only.
        # NOT used as a filter.
        ####################################################

        global_alt_af = (
            total_alt
            / total_an
        )


        ####################################################
        # Retain
        ####################################################

        mrk += 1

        qc[
            "retained_polymorphic_SVs"
        ] += 1


        ####################################################
        # BayPass geno:
        #
        # cx_REF cx_ALT
        # fj_REF fj_ALT
        # hb_REF hb_ALT
        # hn_REF hn_ALT
        # zm_REF zm_ALT
        ####################################################

        geno_values = []


        for p in POP_ORDER:

            geno_values.extend([
                str(ref_count[p]),
                str(alt_count[p])
            ])


        geno_out.write(
            " ".join(geno_values)
            + "\n"
        )


        ####################################################
        # Marker metadata
        ####################################################

        row = [
            str(mrk),
            chrom,
            pos,
            vid,
            ref,
            alt,
            filt,
            str(end),
            str(svtype),
            str(svlen),
            str(total_ref),
            str(total_alt),
            str(total_an),
            f"{global_alt_af:.8f}"
        ]


        for p in POP_ORDER:

            an = (
                ref_count[p]
                + alt_count[p]
            )


            if an > 0:

                af = (
                    alt_count[p]
                    / an
                )

            else:

                af = float("nan")


            row.extend([
                str(ref_count[p]),
                str(alt_count[p]),
                str(an),
                f"{af:.8f}",
                str(n_called[p])
            ])


        marker_out.write(
            "\t".join(row)
            + "\n"
        )


############################################################
# Close
############################################################

geno_out.close()
marker_out.close()


############################################################
# Invalid genotype report
############################################################

with open(
    INVALID_GT_OUT,
    "w"
) as out:

    out.write(
        "CHROM\tPOS\tID\tSAMPLE\tRAW_GT\n"
    )

    for row in invalid_gt_records:

        out.write(
            "\t".join(
                map(str, row)
            )
            + "\n"
        )


############################################################
# QC report
############################################################

with open(
    QC_OUT,
    "w"
) as out:

    out.write(
        "metric\tcount\n"
    )

    for key, value in qc.items():

        out.write(
            f"{key}\t{value}\n"
        )


############################################################
# Console summary
############################################################

print("\nQC summary")

for key, value in qc.items():

    print(
        key,
        value,
        sep="\t"
    )


print(
    "\nOutput geno:",
    GENO_OUT
)

print(
    "Output marker:",
    MARKER_OUT
)

print(
    "Invalid GT report:",
    INVALID_GT_OUT
)

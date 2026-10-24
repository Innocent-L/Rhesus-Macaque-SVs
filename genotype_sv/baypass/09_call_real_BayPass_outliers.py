#!/usr/bin/env python3

import csv
import sys
import math


############################################################
# Usage
############################################################

if len(sys.argv) != 3:

    sys.exit(
        "\nUsage:\n"
        "python3 09_call_real_BayPass_outliers.py "
        "china46_core_XtX_with_SV_ID.tsv "
        "04_pod/POD_XtX_thresholds.tsv\n"
    )


REAL_FILE = sys.argv[1]
THRESHOLD_FILE = sys.argv[2]


############################################################
# 1. Read POD thresholds
############################################################

thresholds = {}


with open(THRESHOLD_FILE, "r") as f:

    for line_no, line in enumerate(f, start=1):

        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 2:

            raise RuntimeError(
                f"ERROR: malformed threshold file "
                f"at line {line_no}: {line.strip()}"
            )

        key, value = fields

        thresholds[key] = value


required_thresholds = [
    "POD95_M_XtX",
    "POD975_M_XtX",
    "POD99_M_XtX"
]


for key in required_thresholds:

    if key not in thresholds:

        raise RuntimeError(
            f"ERROR: required threshold "
            f"'{key}' was not found."
        )


T95 = float(
    thresholds["POD95_M_XtX"]
)

T975 = float(
    thresholds["POD975_M_XtX"]
)

T99 = float(
    thresholds["POD99_M_XtX"]
)


############################################################
# Check threshold ordering
############################################################

if not (
    T95 <= T975 <= T99
):
    raise RuntimeError(
        "ERROR: threshold ordering is invalid: "
        f"POD95={T95}, "
        f"POD97.5={T975}, "
        f"POD99={T99}"
    )


print(
    "POD95 threshold:",
    T95
)

print(
    "POD97.5 threshold:",
    T975
)

print(
    "POD99 threshold:",
    T99
)


############################################################
# 2. Read real BayPass results
############################################################

with open(
    REAL_FILE,
    "r",
    newline=""
) as f:

    reader = csv.DictReader(
        f,
        delimiter="\t"
    )

    fields = reader.fieldnames

    if fields is None:
        raise RuntimeError(
            "ERROR: no header found "
            "in the real BayPass result file."
        )


    required_columns = [
        "M_XtX",
        "CHROM",
        "POS"
    ]


    for col in required_columns:

        if col not in fields:

            raise RuntimeError(
                f"ERROR: required column "
                f"'{col}' was not found "
                f"in {REAL_FILE}"
            )


    rows = list(reader)


############################################################
# 3. Mark outlier levels
############################################################

valid_rows = []


for i, row in enumerate(
    rows,
    start=2
):

    try:

        xtx = float(
            row["M_XtX"]
        )

    except ValueError:

        raise RuntimeError(
            f"ERROR: invalid M_XtX "
            f"at line {i}: "
            f"{row['M_XtX']}"
        )


    if not math.isfinite(xtx):

        raise RuntimeError(
            f"ERROR: non-finite M_XtX "
            f"at line {i}"
        )


    ########################################################
    # Threshold columns
    ########################################################

    row[
        "POD95_threshold"
    ] = f"{T95:.10f}"

    row[
        "POD975_threshold"
    ] = f"{T975:.10f}"

    row[
        "POD99_threshold"
    ] = f"{T99:.10f}"


    ########################################################
    # Outlier flags
    #
    # >= is used so loci exactly equal to a POD
    # quantile are retained.
    ########################################################

    row[
        "BayPass_POD95"
    ] = (
        "TRUE"
        if xtx >= T95
        else "FALSE"
    )


    row[
        "BayPass_POD975"
    ] = (
        "TRUE"
        if xtx >= T975
        else "FALSE"
    )


    row[
        "BayPass_POD99"
    ] = (
        "TRUE"
        if xtx >= T99
        else "FALSE"
    )


    ########################################################
    # Highest support level reached
    #
    # This gives mutually exclusive categories:
    #
    # POD99
    # POD97.5_only
    # POD95_only
    # Below_POD95
    ########################################################

    if xtx >= T99:

        row[
            "BayPass_support_level"
        ] = "POD99"

    elif xtx >= T975:

        row[
            "BayPass_support_level"
        ] = "POD97.5_only"

    elif xtx >= T95:

        row[
            "BayPass_support_level"
        ] = "POD95_only"

    else:

        row[
            "BayPass_support_level"
        ] = "Below_POD95"


    valid_rows.append(
        row
    )


############################################################
# 4. Output field names
############################################################

new_fields = (
    fields
    + [
        "POD95_threshold",
        "POD975_threshold",
        "POD99_threshold",
        "BayPass_POD95",
        "BayPass_POD975",
        "BayPass_POD99",
        "BayPass_support_level"
    ]
)


############################################################
# 5. Full calibrated result
############################################################

FULL_OUT = (
    "china46_core_XtX_PODcalibrated.tsv"
)


with open(
    FULL_OUT,
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=new_fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        valid_rows
    )


############################################################
# 6. Define nested outlier sets
############################################################

out95 = [
    r
    for r in valid_rows
    if r["BayPass_POD95"]
    == "TRUE"
]


out975 = [
    r
    for r in valid_rows
    if r["BayPass_POD975"]
    == "TRUE"
]


out99 = [
    r
    for r in valid_rows
    if r["BayPass_POD99"]
    == "TRUE"
]


############################################################
# Sort by decreasing XtX
############################################################

def sort_by_xtx(rows):

    return sorted(
        rows,
        key=lambda r:
        float(
            r["M_XtX"]
        ),
        reverse=True
    )


out95 = sort_by_xtx(
    out95
)

out975 = sort_by_xtx(
    out975
)

out99 = sort_by_xtx(
    out99
)


############################################################
# 7. Write function
############################################################

def write_table(
    filename,
    data
):

    with open(
        filename,
        "w",
        newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=new_fields,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writeheader()
        writer.writerows(data)


############################################################
# 8. Write three nested outlier datasets
############################################################

write_table(
    "china46_core_POD95_outliers.tsv",
    out95
)


write_table(
    "china46_core_POD975_outliers.tsv",
    out975
)


write_table(
    "china46_core_POD99_outliers.tsv",
    out99
)


############################################################
# 9. Mutually exclusive levels
#
# Useful for reporting:
#
# POD95_only
# POD97.5_only
# POD99
############################################################

pod95_only = [
    r
    for r in valid_rows
    if r[
        "BayPass_support_level"
    ] == "POD95_only"
]


pod975_only = [
    r
    for r in valid_rows
    if r[
        "BayPass_support_level"
    ] == "POD97.5_only"
]


pod99_strict = [
    r
    for r in valid_rows
    if r[
        "BayPass_support_level"
    ] == "POD99"
]


write_table(
    "china46_core_POD95_only.tsv",
    sort_by_xtx(
        pod95_only
    )
)


write_table(
    "china46_core_POD975_only.tsv",
    sort_by_xtx(
        pod975_only
    )
)


write_table(
    "china46_core_POD99_strict.tsv",
    sort_by_xtx(
        pod99_strict
    )
)


############################################################
# 10. Check nesting
#
# Use CHROM+POS rather than ID,
# because your SV IDs are not unique.
############################################################

def locus_key(row):

    return (
        row["CHROM"],
        row["POS"]
    )


set95 = {
    locus_key(r)
    for r in out95
}


set975 = {
    locus_key(r)
    for r in out975
}


set99 = {
    locus_key(r)
    for r in out99
}


if not set99.issubset(
    set975
):

    raise RuntimeError(
        "ERROR: POD99 set is not "
        "a subset of POD97.5."
    )


if not set975.issubset(
    set95
):

    raise RuntimeError(
        "ERROR: POD97.5 set is not "
        "a subset of POD95."
    )


############################################################
# 11. Summary statistics
############################################################

N = len(
    valid_rows
)

N95 = len(
    out95
)

N975 = len(
    out975
)

N99 = len(
    out99
)


############################################################
# Mutually exclusive numbers
############################################################

N95_only = len(
    pod95_only
)

N975_only = len(
    pod975_only
)

N99_strict = len(
    pod99_strict
)


############################################################
# 12. Write summary table
############################################################

SUMMARY_OUT = (
    "china46_core_POD_outlier_summary.tsv"
)


with open(
    SUMMARY_OUT,
    "w"
) as f:

    f.write(
        "metric\tcount\tproportion\n"
    )


    f.write(
        f"Total_real_SVs\t"
        f"{N}\t"
        f"1.000000\n"
    )


    f.write(
        f"POD95_outliers\t"
        f"{N95}\t"
        f"{N95/N:.6f}\n"
    )


    f.write(
        f"POD97.5_outliers\t"
        f"{N975}\t"
        f"{N975/N:.6f}\n"
    )


    f.write(
        f"POD99_outliers\t"
        f"{N99}\t"
        f"{N99/N:.6f}\n"
    )


    f.write(
        f"POD95_only\t"
        f"{N95_only}\t"
        f"{N95_only/N:.6f}\n"
    )


    f.write(
        f"POD97.5_only\t"
        f"{N975_only}\t"
        f"{N975_only/N:.6f}\n"
    )


    f.write(
        f"POD99_strict\t"
        f"{N99_strict}\t"
        f"{N99_strict/N:.6f}\n"
    )


############################################################
# 13. Console summary
############################################################

print()
print(
    "========================================"
)

print(
    "BayPass POD-calibrated XtX summary"
)

print(
    "========================================"
)


print(
    "Total real SVs:",
    N
)


print(
    "POD95 threshold:",
    f"{T95:.10f}"
)

print(
    "POD95 outliers:",
    N95,
    f"({N95/N:.4%})"
)


print(
    "POD97.5 threshold:",
    f"{T975:.10f}"
)

print(
    "POD97.5 outliers:",
    N975,
    f"({N975/N:.4%})"
)


print(
    "POD99 threshold:",
    f"{T99:.10f}"
)

print(
    "POD99 outliers:",
    N99,
    f"({N99/N:.4%})"
)


print()
print(
    "Mutually exclusive support levels:"
)

print(
    "POD95 only:",
    N95_only
)

print(
    "POD97.5 only:",
    N975_only
)

print(
    "POD99:",
    N99_strict
)


print()
print(
    "Nested-set QC:"
)

print(
    "POD99 subset POD97.5: TRUE"
)

print(
    "POD97.5 subset POD95: TRUE"
)


print()
print(
    "Output files:"
)

print(
    FULL_OUT
)

print(
    "china46_core_POD95_outliers.tsv"
)

print(
    "china46_core_POD975_outliers.tsv"
)

print(
    "china46_core_POD99_outliers.tsv"
)

print(
    SUMMARY_OUT
)

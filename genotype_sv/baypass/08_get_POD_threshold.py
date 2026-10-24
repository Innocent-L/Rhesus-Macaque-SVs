#!/usr/bin/env python3

import sys
import math


############################################################
# Usage
############################################################

if len(sys.argv) != 2:
    sys.exit(
        "\nUsage:\n"
        "python3 08_get_POD_threshold.py "
        "04_pod/china46_POD_summary_pi_xtx.out\n"
    )


INPUT = sys.argv[1]


############################################################
# Read POD M_XtX values
############################################################

values = []

with open(INPUT, "r") as f:

    header = f.readline().split()

    if "M_XtX" not in header:
        raise RuntimeError(
            "ERROR: column 'M_XtX' was not found "
            "in the POD summary file."
        )

    xtx_idx = header.index("M_XtX")


    for line_no, line in enumerate(f, start=2):

        if not line.strip():
            continue

        fields = line.split()

        if len(fields) <= xtx_idx:
            raise RuntimeError(
                f"ERROR: malformed line {line_no}: "
                f"not enough columns."
            )

        try:
            value = float(
                fields[xtx_idx]
            )

        except ValueError:
            raise RuntimeError(
                f"ERROR: invalid M_XtX value "
                f"at line {line_no}: "
                f"{fields[xtx_idx]}"
            )


        if math.isfinite(value):
            values.append(value)


############################################################
# Basic QC
############################################################

if len(values) == 0:
    raise RuntimeError(
        "ERROR: no valid M_XtX values were read."
    )


values.sort()

N = len(values)


############################################################
# R quantile(type = 7) equivalent
#
# R default:
# quantile(x, probs, type = 7)
############################################################

def quantile_type7(data, p):

    n = len(data)

    if n == 0:
        raise ValueError(
            "Cannot calculate quantile of empty data."
        )

    if n == 1:
        return data[0]


    h = (n - 1) * p

    lower = math.floor(h)
    upper = math.ceil(h)

    fraction = h - lower


    if lower == upper:
        return data[lower]


    return (
        data[lower]
        +
        fraction
        *
        (
            data[upper]
            -
            data[lower]
        )
    )


############################################################
# POD thresholds
############################################################

q95 = quantile_type7(
    values,
    0.95
)

q975 = quantile_type7(
    values,
    0.975
)

q99 = quantile_type7(
    values,
    0.99
)


############################################################
# Output
############################################################

print(
    f"POD_loci\t{N}"
)

print(
    f"POD95_M_XtX\t{q95:.10f}"
)

print(
    f"POD975_M_XtX\t{q975:.10f}"
)

print(
    f"POD99_M_XtX\t{q99:.10f}"
)

print(
    f"POD_min_M_XtX\t{values[0]:.10f}"
)

print(
    f"POD_max_M_XtX\t{values[-1]:.10f}"
)

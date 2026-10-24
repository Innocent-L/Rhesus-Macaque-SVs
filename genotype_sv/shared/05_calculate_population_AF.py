import pandas as pd
import numpy as np

gt_file = "shared_population_noLR_GT.tsv"
map_file = "sample_population.tsv"

# ============================================================
# Read data
# ============================================================

gt = pd.read_csv(
    gt_file,
    sep="\t",
    dtype=str
)

mapping = pd.read_csv(
    map_file,
    sep="\t",
    header=None,
    names=["sample", "population"],
    dtype=str
)

meta_cols = [
    "CHROM",
    "POS",
    "END",
    "ID",
    "SVTYPE",
    "SVLEN"
]

sample_cols = [
    c for c in gt.columns
    if c not in meta_cols
]

print("Number of SVs:", len(gt))
print("Number of samples:", len(sample_cols))

# Check mapping
missing_mapping = set(sample_cols) - set(mapping["sample"])

if missing_mapping:
    print("ERROR: samples missing from population mapping:")
    for x in sorted(missing_mapping):
        print(x)
    raise SystemExit

mapping = mapping[
    mapping["sample"].isin(sample_cols)
].copy()

print("\nSamples per population:")
print(mapping["population"].value_counts().sort_index())


# ============================================================
# Genotype parser
# ============================================================

def parse_gt(x):
    """
    Returns:
        AC contribution
        AN contribution
        carrier status

    Examples:
        0/0 -> AC=0, AN=2, carrier=0
        0/1 -> AC=1, AN=2, carrier=1
        1/1 -> AC=2, AN=2, carrier=1
        ./. -> AC=0, AN=0, carrier=NA
    """

    if pd.isna(x):
        return 0, 0, np.nan

    x = str(x)

    if ":" in x:
        x = x.split(":")[0]

    if x in [".", "./.", ".|.", ""]:
        return 0, 0, np.nan

    alleles = x.replace("|", "/").split("/")

    called = [
        a for a in alleles
        if a != "."
    ]

    if len(called) == 0:
        return 0, 0, np.nan

    ac = 0

    for a in called:
        try:
            allele = int(a)

            if allele > 0:
                ac += 1

        except ValueError:
            pass

    an = len(called)

    carrier = 1 if ac > 0 else 0

    return ac, an, carrier


# ============================================================
# Calculate population-specific AC / AN / AF
# ============================================================

results = []

populations = sorted(
    mapping["population"].unique()
)

for i, row in gt.iterrows():

    base = {
        "CHROM": row["CHROM"],
        "POS": row["POS"],
        "END": row["END"],
        "ID": row["ID"],
        "SVTYPE": row["SVTYPE"],
        "SVLEN": row["SVLEN"]
    }

    for pop in populations:

        samples = mapping.loc[
            mapping["population"] == pop,
            "sample"
        ].tolist()

        total_ac = 0
        total_an = 0
        n_called = 0
        n_carriers = 0

        for sample in samples:

            ac, an, carrier = parse_gt(
                row[sample]
            )

            total_ac += ac
            total_an += an

            if an > 0:
                n_called += 1

            if carrier == 1:
                n_carriers += 1

        if total_an > 0:
            af = total_ac / total_an
        else:
            af = np.nan

        results.append({
            **base,
            "Population": pop,
            "N_samples_total": len(samples),
            "N_samples_called": n_called,
            "N_carriers": n_carriers,
            "AC": total_ac,
            "AN": total_an,
            "AF": af,
            "Present": (
                1 if total_ac > 0
                else 0 if total_an > 0
                else np.nan
            )
        })


long = pd.DataFrame(results)

long.to_csv(
    "shared_SV_population_AC_AN_AF_long.tsv",
    sep="\t",
    index=False
)

print("\nLong table written.")
print("Rows:", len(long))


# ============================================================
# Population summary
# ============================================================

summary = []

for pop in populations:

    x = long[
        long["Population"] == pop
    ].copy()

    callable_sites = x[
        x["AN"] > 0
    ]

    present_sites = callable_sites[
        callable_sites["AC"] > 0
    ]

    summary.append({
        "Population": pop,

        "Total_shared_SVs": len(x),

        "Callable_SVs":
            len(callable_sites),

        "Present_SVs_AC_gt0":
            len(present_sites),

        "Absent_SVs_AC_eq0":
            len(callable_sites) -
            len(present_sites),

        "Proportion_present_among_callable":
            (
                len(present_sites) /
                len(callable_sites)
                if len(callable_sites) > 0
                else np.nan
            ),

        # IMPORTANT:
        # zeros are retained here
        "Mean_AF_all_callable":
            callable_sites["AF"].mean(),

        "Median_AF_all_callable":
            callable_sites["AF"].median(),

        # descriptive only among present variants
        "Mean_AF_present_only":
            present_sites["AF"].mean()
            if len(present_sites) > 0
            else np.nan,

        "Median_AF_present_only":
            present_sites["AF"].median()
            if len(present_sites) > 0
            else np.nan,

        "Mean_called_samples_per_SV":
            callable_sites["N_samples_called"].mean()
    })


summary = pd.DataFrame(summary)

summary.to_csv(
    "shared_SV_population_summary.tsv",
    sep="\t",
    index=False
)

print("\nPopulation summary:")
print(summary.to_string(index=False))


# ============================================================
# Wide AF matrix
# One row = one SV
# One column = one population
# ============================================================

long["SV_KEY"] = (
    long["CHROM"].astype(str) + ":" +
    long["POS"].astype(str) + ":" +
    long["END"].astype(str) + ":" +
    long["SVTYPE"].astype(str) + ":" +
    long["SVLEN"].astype(str)
)

af_wide = long.pivot_table(
    index=[
        "SV_KEY",
        "CHROM",
        "POS",
        "END",
        "ID",
        "SVTYPE",
        "SVLEN"
    ],
    columns="Population",
    values="AF",
    aggfunc="first"
).reset_index()

af_wide.to_csv(
    "shared_SV_population_AF_wide.tsv",
    sep="\t",
    index=False
)


# ============================================================
# Presence/absence wide matrix
# ============================================================

presence_wide = long.pivot_table(
    index=[
        "SV_KEY",
        "CHROM",
        "POS",
        "END",
        "ID",
        "SVTYPE",
        "SVLEN"
    ],
    columns="Population",
    values="Present",
    aggfunc="first"
).reset_index()

presence_wide.to_csv(
    "shared_SV_population_presence_wide.tsv",
    sep="\t",
    index=False
)


print("\nFinished.")

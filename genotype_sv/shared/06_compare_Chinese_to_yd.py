import pandas as pd
import numpy as np

# ============================================================
# Settings
# ============================================================

infile = "shared_SV_population_AC_AN_AF_long.tsv"

# Indian-origin population
indian_pop = "yd"

# Five Chinese populations
chinese_pops = [
    "cx",
    "fj",
    "hb",
    "hn",
    "zm"
]

all_pops = chinese_pops + [indian_pop]

# Bootstrap settings
N_BOOT = 5000
SEED = 20260926


# ============================================================
# Spearman correlation without scipy
# ============================================================

def spearman_corr(x, y):
    """
    Calculate Spearman rank correlation without scipy.

    Spearman correlation is equivalent to Pearson correlation
    calculated on ranked values.

    Ties are handled using average ranks.
    Missing values are removed pairwise.
    """

    x = pd.Series(x, dtype=float)
    y = pd.Series(y, dtype=float)

    valid = x.notna() & y.notna()

    x = x[valid]
    y = y[valid]

    if len(x) < 2:
        return np.nan

    rx = x.rank(method="average")
    ry = y.rank(method="average")

    # Avoid undefined correlation for constant vectors
    if rx.nunique() <= 1 or ry.nunique() <= 1:
        return np.nan

    return np.corrcoef(
        rx.to_numpy(),
        ry.to_numpy()
    )[0, 1]


# ============================================================
# Read input
# ============================================================

df = pd.read_csv(
    infile,
    sep="\t"
)

print("================================")
print("INPUT")
print("================================")

print("Input rows:", len(df))

print(
    "Unique populations:",
    sorted(df["Population"].unique())
)


# ============================================================
# Check required columns
# ============================================================

required_cols = [
    "CHROM",
    "POS",
    "END",
    "ID",
    "SVTYPE",
    "SVLEN",
    "Population",
    "AC",
    "AN",
    "AF"
]

missing_cols = [
    c for c in required_cols
    if c not in df.columns
]

if missing_cols:
    raise ValueError(
        f"Missing required columns: {missing_cols}"
    )


# ============================================================
# Check populations
# ============================================================

available_pops = set(
    df["Population"].unique()
)

missing_pops = [
    p for p in all_pops
    if p not in available_pops
]

if missing_pops:
    raise ValueError(
        f"Required populations missing from input: {missing_pops}"
    )


# ============================================================
# Build stable SV key
# ============================================================

df["SV_KEY"] = (
    df["CHROM"].astype(str) + ":" +
    df["POS"].astype(str) + ":" +
    df["END"].astype(str) + ":" +
    df["SVTYPE"].astype(str) + ":" +
    df["SVLEN"].astype(str)
)


print(
    "Unique SVs:",
    df["SV_KEY"].nunique()
)


# ============================================================
# Keep only populations of interest
# ============================================================

x = df[
    df["Population"].isin(all_pops)
].copy()


# ============================================================
# Convert numeric columns safely
# ============================================================

for col in ["AC", "AN", "AF"]:
    x[col] = pd.to_numeric(
        x[col],
        errors="coerce"
    )


# ============================================================
# Identify common-callable SV set
#
# Definition:
# AN > 0 in all six populations
# ============================================================

callable_matrix = x.pivot_table(
    index="SV_KEY",
    columns="Population",
    values="AN",
    aggfunc="first"
)

# Check all six populations are present as columns
for pop in all_pops:
    if pop not in callable_matrix.columns:
        raise ValueError(
            f"{pop} is absent from callable matrix"
        )

common_mask = (
    callable_matrix[all_pops] > 0
).all(axis=1)

common_keys = callable_matrix.index[
    common_mask
].tolist()


total_shared = df["SV_KEY"].nunique()
n_common = len(common_keys)

print()
print("================================")
print("COMMON-CALLABLE SET")
print("================================")

print(
    "Total population-QC shared SVs:",
    total_shared
)

print(
    "Common-callable SVs:",
    n_common
)

print(
    "Retention proportion:",
    n_common / total_shared
)


# Save common-callable keys
pd.Series(
    common_keys,
    name="SV_KEY"
).to_csv(
    "common_callable_shared_SVKEY.txt",
    index=False,
    header=False
)


# ============================================================
# Restrict to common-callable SVs
# ============================================================

common = x[
    x["SV_KEY"].isin(common_keys)
].copy()

common.to_csv(
    "shared_SV_common_callable_long.tsv",
    sep="\t",
    index=False
)


# ============================================================
# Generate AF matrix
# ============================================================

af = common.pivot_table(
    index="SV_KEY",
    columns="Population",
    values="AF",
    aggfunc="first"
)

af = af[all_pops]


# ============================================================
# Generate AC matrix
# ============================================================

ac = common.pivot_table(
    index="SV_KEY",
    columns="Population",
    values="AC",
    aggfunc="first"
)

ac = ac[all_pops]


# ============================================================
# Check AF missingness
# ============================================================

print()
print("AF missing values in common-callable set:")

for pop in all_pops:
    print(
        pop,
        int(af[pop].isna().sum())
    )


# A callable site should normally have AF.
# Keep only sites with non-missing AF in all populations.
af_complete_mask = (
    af[all_pops]
    .notna()
    .all(axis=1)
)

complete_keys = af.index[
    af_complete_mask
]

af = af.loc[
    complete_keys
].copy()

ac = ac.loc[
    complete_keys
].copy()


print()
print(
    "Common-callable SVs with complete AF:",
    len(af)
)


# ============================================================
# Presence / absence
#
# Present = AC > 0
# ============================================================

presence = (
    ac > 0
).astype(int)


# ============================================================
# Save AF and presence matrices
# ============================================================

af.to_csv(
    "common_callable_AF_matrix.tsv",
    sep="\t"
)

presence.to_csv(
    "common_callable_presence_matrix.tsv",
    sep="\t"
)


# ============================================================
# Population-level summary
# ============================================================

summary_rows = []

for pop in all_pops:

    n_present = int(
        presence[pop].sum()
    )

    n_absent = len(presence) - n_present

    summary_rows.append({
        "Population": pop,
        "Common_callable_SVs": len(presence),
        "Present_SVs": n_present,
        "Absent_SVs": n_absent,
        "Proportion_present":
            n_present / len(presence),
        "Mean_AF":
            af[pop].mean(),
        "Median_AF":
            af[pop].median()
    })


summary = pd.DataFrame(
    summary_rows
)

summary.to_csv(
    "common_callable_population_summary.tsv",
    sep="\t",
    index=False
)


print()
print("================================")
print("COMMON-CALLABLE POPULATION SUMMARY")
print("================================")

print(
    summary.to_string(index=False)
)


# ============================================================
# Compare each Chinese population with yd
# ============================================================

yd_af = af[indian_pop].to_numpy(
    dtype=float
)

yd_presence = presence[indian_pop].to_numpy(
    dtype=int
)


similarity_rows = []


for pop in chinese_pops:

    pop_af = af[pop].to_numpy(
        dtype=float
    )

    pop_presence = presence[pop].to_numpy(
        dtype=int
    )

    # --------------------------------------------------------
    # Allele-frequency difference
    # --------------------------------------------------------

    abs_diff = np.abs(
        pop_af - yd_af
    )

    mean_abs_diff = np.mean(
        abs_diff
    )

    median_abs_diff = np.median(
        abs_diff
    )

    # --------------------------------------------------------
    # Root mean squared error
    # --------------------------------------------------------

    rmse = np.sqrt(
        np.mean(
            (pop_af - yd_af) ** 2
        )
    )

    # --------------------------------------------------------
    # Spearman rank correlation
    # --------------------------------------------------------

    rho = spearman_corr(
        pop_af,
        yd_af
    )

    # --------------------------------------------------------
    # Presence / absence statistics
    # --------------------------------------------------------

    both_present = np.sum(
        (pop_presence == 1) &
        (yd_presence == 1)
    )

    china_only = np.sum(
        (pop_presence == 1) &
        (yd_presence == 0)
    )

    yd_only = np.sum(
        (pop_presence == 0) &
        (yd_presence == 1)
    )

    both_absent = np.sum(
        (pop_presence == 0) &
        (yd_presence == 0)
    )

    union_present = np.sum(
        (pop_presence == 1) |
        (yd_presence == 1)
    )

    if union_present > 0:
        jaccard = (
            both_present /
            union_present
        )
    else:
        jaccard = np.nan


    similarity_rows.append({
        "Chinese_population": pop,
        "N_common_callable": len(af),

        "China_present":
            int(pop_presence.sum()),

        "YD_present":
            int(yd_presence.sum()),

        "Both_present":
            int(both_present),

        "China_only":
            int(china_only),

        "YD_only":
            int(yd_only),

        "Both_absent":
            int(both_absent),

        "Jaccard_presence":
            jaccard,

        "Mean_abs_AF_difference":
            mean_abs_diff,

        "Median_abs_AF_difference":
            median_abs_diff,

        "RMSE_AF":
            rmse,

        "Spearman_rho":
            rho
    })


similarity = pd.DataFrame(
    similarity_rows
)


# ============================================================
# Bootstrap confidence intervals
#
# Resample SVs with replacement.
# ============================================================

rng = np.random.default_rng(
    SEED
)

n = len(af)

bootstrap_CI = {}


print()
print("================================")
print("BOOTSTRAP")
print("================================")

print(
    "Bootstrap iterations:",
    N_BOOT
)

print(
    "Number of SVs per bootstrap:",
    n
)


for pop in chinese_pops:

    print(
        f"Running bootstrap: {pop} vs {indian_pop}"
    )

    pop_af = af[pop].to_numpy(
        dtype=float
    )

    pop_presence = presence[pop].to_numpy(
        dtype=int
    )

    mad_boot = np.empty(
        N_BOOT,
        dtype=float
    )

    rho_boot = np.empty(
        N_BOOT,
        dtype=float
    )

    jaccard_boot = np.empty(
        N_BOOT,
        dtype=float
    )


    for b in range(N_BOOT):

        idx = rng.integers(
            low=0,
            high=n,
            size=n
        )

        xaf = pop_af[idx]
        yaf = yd_af[idx]

        xp = pop_presence[idx]
        yp = yd_presence[idx]


        # Mean absolute AF difference
        mad_boot[b] = np.mean(
            np.abs(
                xaf - yaf
            )
        )


        # Spearman correlation
        rho_boot[b] = spearman_corr(
            xaf,
            yaf
        )


        # Jaccard
        intersection = np.sum(
            (xp == 1) &
            (yp == 1)
        )

        union = np.sum(
            (xp == 1) |
            (yp == 1)
        )

        if union > 0:
            jaccard_boot[b] = (
                intersection /
                union
            )
        else:
            jaccard_boot[b] = np.nan


    bootstrap_CI[pop] = {

        "MAD_low":
            np.nanpercentile(
                mad_boot,
                2.5
            ),

        "MAD_high":
            np.nanpercentile(
                mad_boot,
                97.5
            ),

        "Rho_low":
            np.nanpercentile(
                rho_boot,
                2.5
            ),

        "Rho_high":
            np.nanpercentile(
                rho_boot,
                97.5
            ),

        "Jaccard_low":
            np.nanpercentile(
                jaccard_boot,
                2.5
            ),

        "Jaccard_high":
            np.nanpercentile(
                jaccard_boot,
                97.5
            )
    }


# ============================================================
# Add bootstrap CIs to final table
# ============================================================

similarity[
    "Mean_abs_AF_difference_CI_low"
] = np.nan

similarity[
    "Mean_abs_AF_difference_CI_high"
] = np.nan

similarity[
    "Spearman_rho_CI_low"
] = np.nan

similarity[
    "Spearman_rho_CI_high"
] = np.nan

similarity[
    "Jaccard_CI_low"
] = np.nan

similarity[
    "Jaccard_CI_high"
] = np.nan


for pop in chinese_pops:

    mask = (
        similarity["Chinese_population"]
        == pop
    )

    similarity.loc[
        mask,
        "Mean_abs_AF_difference_CI_low"
    ] = bootstrap_CI[pop]["MAD_low"]

    similarity.loc[
        mask,
        "Mean_abs_AF_difference_CI_high"
    ] = bootstrap_CI[pop]["MAD_high"]

    similarity.loc[
        mask,
        "Spearman_rho_CI_low"
    ] = bootstrap_CI[pop]["Rho_low"]

    similarity.loc[
        mask,
        "Spearman_rho_CI_high"
    ] = bootstrap_CI[pop]["Rho_high"]

    similarity.loc[
        mask,
        "Jaccard_CI_low"
    ] = bootstrap_CI[pop]["Jaccard_low"]

    similarity.loc[
        mask,
        "Jaccard_CI_high"
    ] = bootstrap_CI[pop]["Jaccard_high"]


# ============================================================
# Save final similarity table
# ============================================================

similarity.to_csv(
    "Chinese_vs_yd_shared_SV_similarity.tsv",
    sep="\t",
    index=False
)


# ============================================================
# Also save a compact reviewer-friendly table
# ============================================================

compact_cols = [
    "Chinese_population",
    "N_common_callable",
    "Both_present",
    "Jaccard_presence",
    "Jaccard_CI_low",
    "Jaccard_CI_high",
    "Mean_abs_AF_difference",
    "Mean_abs_AF_difference_CI_low",
    "Mean_abs_AF_difference_CI_high",
    "Median_abs_AF_difference",
    "RMSE_AF",
    "Spearman_rho",
    "Spearman_rho_CI_low",
    "Spearman_rho_CI_high"
]

compact = similarity[
    compact_cols
].copy()

compact.to_csv(
    "Chinese_vs_yd_shared_SV_similarity_compact.tsv",
    sep="\t",
    index=False
)


# ============================================================
# Print final result
# ============================================================

print()
print("================================")
print("CHINESE vs YD SIMILARITY")
print("================================")

print(
    compact.to_string(
        index=False
    )
)

print()
print("Finished successfully.")

print()
print("Main output files:")
print(
    "1. common_callable_population_summary.tsv"
)
print(
    "2. Chinese_vs_yd_shared_SV_similarity.tsv"
)
print(
    "3. Chinese_vs_yd_shared_SV_similarity_compact.tsv"
)
print(
    "4. common_callable_AF_matrix.tsv"
)
print(
    "5. common_callable_presence_matrix.tsv"
)

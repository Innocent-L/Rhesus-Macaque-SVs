import pandas as pd
import numpy as np


shared_file = "India_China_shared_assembly_full.tsv"
pop_file = "population_QC_all_sites.tsv"
afgt0_file = "population_AFgt0_sites.tsv"


shared = pd.read_csv(shared_file, sep="\t", dtype=str)
pop = pd.read_csv(pop_file, sep="\t", dtype=str)
afgt0 = pd.read_csv(afgt0_file, sep="\t", dtype=str)


def make_key5(df):
    return (
        df["CHROM"].astype(str) + ":" +
        df["POS"].astype(str) + ":" +
        df["END"].astype(str) + ":" +
        df["SVTYPE"].astype(str) + ":" +
        df["SVLEN"].astype(str)
    )


def make_key7(df):
    return (
        make_key5(df) + ":" +
        df["REF"].astype(str) + ":" +
        df["ALT"].astype(str)
    )


for x in [shared, pop, afgt0]:
    x["KEY5"] = make_key5(x)
    x["KEY7"] = make_key7(x)


print("========== INPUT ==========")

print("Assembly shared:", len(shared))
print("Shared unique KEY5:", shared["KEY5"].nunique())
print("Shared unique KEY7:", shared["KEY7"].nunique())

print()

print("Population QC records:", len(pop))
print("Population QC unique KEY5:", pop["KEY5"].nunique())
print("Population QC unique KEY7:", pop["KEY7"].nunique())

print()

print("AF>0 records:", len(afgt0))
print("AF>0 unique KEY5:", afgt0["KEY5"].nunique())
print("AF>0 unique KEY7:", afgt0["KEY7"].nunique())


# ============================================================
# Matching function
#
# Priority:
# 1. Unique KEY5 in both datasets
# 2. For duplicated KEY5, exact KEY7
# ============================================================

def match_datasets(source, target, target_name):

    source = source.copy()
    target = target.copy()

    source_key5_counts = source["KEY5"].value_counts()
    target_key5_counts = target["KEY5"].value_counts()

    target_key7_counts = target["KEY7"].value_counts()

    # Fast lookup
    target_by_key5 = {
        k: g.iloc[0]
        for k, g in target.groupby("KEY5")
        if len(g) == 1
    }

    target_by_key7 = {
        k: g.iloc[0]
        for k, g in target.groupby("KEY7")
        if len(g) == 1
    }

    results = []

    for idx, row in source.iterrows():

        k5 = row["KEY5"]
        k7 = row["KEY7"]

        matched_row = None
        method = "unmatched"

        # ----------------------------------------------------
        # Method 1:
        # KEY5 must be unique in BOTH source and target
        # ----------------------------------------------------
        if (
            source_key5_counts.get(k5, 0) == 1
            and
            target_key5_counts.get(k5, 0) == 1
        ):
            matched_row = target_by_key5[k5]
            method = "KEY5_unique"

        # ----------------------------------------------------
        # Method 2:
        # resolve duplicated coordinate/length records
        # using REF + ALT
        # ----------------------------------------------------
        elif target_key7_counts.get(k7, 0) == 1:
            matched_row = target_by_key7[k7]
            method = "KEY7_strict"

        out = row.to_dict()
        out["match_method"] = method

        if matched_row is not None:
            out["POP_CHROM"] = matched_row["CHROM"]
            out["POP_POS"] = matched_row["POS"]
            out["POP_END"] = matched_row["END"]
            out["POP_ID"] = matched_row["ID"]
            out["POP_SVTYPE"] = matched_row["SVTYPE"]
            out["POP_SVLEN"] = matched_row["SVLEN"]
            out["POP_REF"] = matched_row["REF"]
            out["POP_ALT"] = matched_row["ALT"]
            out["POP_AF"] = matched_row.get("AF", np.nan)
            out["POP_KEY7"] = matched_row["KEY7"]
        else:
            out["POP_CHROM"] = np.nan
            out["POP_POS"] = np.nan
            out["POP_END"] = np.nan
            out["POP_ID"] = np.nan
            out["POP_SVTYPE"] = np.nan
            out["POP_SVLEN"] = np.nan
            out["POP_REF"] = np.nan
            out["POP_ALT"] = np.nan
            out["POP_AF"] = np.nan
            out["POP_KEY7"] = np.nan

        results.append(out)

    result = pd.DataFrame(results)

    result.to_csv(
        f"shared_vs_{target_name}_all.tsv",
        sep="\t",
        index=False
    )

    matched = result[
        result["match_method"] != "unmatched"
    ].copy()

    unmatched = result[
        result["match_method"] == "unmatched"
    ].copy()

    matched.to_csv(
        f"shared_vs_{target_name}_matched.tsv",
        sep="\t",
        index=False
    )

    unmatched.to_csv(
        f"shared_vs_{target_name}_unmatched.tsv",
        sep="\t",
        index=False
    )

    print()
    print(f"========== {target_name} ==========")
    print("Total shared:", len(result))
    print("Matched:", len(matched))
    print("Unmatched:", len(unmatched))
    print("Match rate:", len(matched) / len(result))

    print()
    print("Match methods:")
    print(result["match_method"].value_counts())

    return result


# ============================================================
# 1. Match to complete population-QC VCF
# ============================================================

qc_result = match_datasets(
    shared,
    pop,
    "population_QC"
)


# ============================================================
# 2. Match to AF > 0 VCF
# ============================================================

af_result = match_datasets(
    shared,
    afgt0,
    "population_AFgt0"
)


# ============================================================
# Parse AF from complete QC dataset
# ============================================================

def parse_af(x):

    if pd.isna(x):
        return np.nan

    x = str(x)

    if x in [".", "", "nan", "NaN"]:
        return np.nan

    try:
        vals = [
            float(i)
            for i in x.split(",")
            if i not in ["", "."]
        ]

        if len(vals) == 0:
            return np.nan

        return max(vals)

    except:
        return np.nan


qc_result["POP_AF_numeric"] = qc_result["POP_AF"].apply(parse_af)


# ============================================================
# Status classification
# ============================================================

qc_result["population_status"] = "not_in_population_QC"

mask_match = qc_result["match_method"] != "unmatched"

qc_result.loc[
    mask_match &
    (qc_result["POP_AF_numeric"] > 0),
    "population_status"
] = "population_QC_AF_gt0"

qc_result.loc[
    mask_match &
    (qc_result["POP_AF_numeric"] == 0),
    "population_status"
] = "population_QC_AF_eq0"

qc_result.loc[
    mask_match &
    qc_result["POP_AF_numeric"].isna(),
    "population_status"
] = "population_QC_AF_missing"


qc_result.to_csv(
    "India_China_shared_population_status.tsv",
    sep="\t",
    index=False
)


status_summary = (
    qc_result["population_status"]
    .value_counts()
    .rename_axis("Status")
    .reset_index(name="N_SVs")
)

status_summary["Proportion_of_21307"] = (
    status_summary["N_SVs"] /
    len(qc_result)
)

status_summary.to_csv(
    "India_China_shared_population_status_summary.tsv",
    sep="\t",
    index=False
)

print()
print("========== FINAL STATUS ==========")
print(status_summary.to_string(index=False))

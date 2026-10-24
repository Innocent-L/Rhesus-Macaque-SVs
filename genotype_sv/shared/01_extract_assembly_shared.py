import pandas as pd

infile = "merge_pan_for_shared.tsv"

df = pd.read_csv(infile, sep="\t", dtype=str)

chinese = ["cx_2", "fj", "hb_1", "xz_1", "zm_2"]
india = "yd_10"


def is_present(gt):
    """
    An SV is considered present if the genotype contains
    any non-reference allele.
    """
    if pd.isna(gt):
        return False

    gt = str(gt).split(":")[0]

    if gt in [".", "./.", ".|."]:
        return False

    alleles = gt.replace("|", "/").split("/")

    return any(a not in ["0", "."] for a in alleles)


# Presence/absence
for sample in chinese + [india]:
    df[sample + "_present"] = (
        df[sample].apply(is_present).astype(int)
    )

# Number of Chinese assemblies carrying each SV
df["N_Chinese"] = df[
    [x + "_present" for x in chinese]
].sum(axis=1)

# India + at least one Chinese assembly
shared = df[
    (df[india + "_present"] == 1) &
    (df["N_Chinese"] >= 1)
].copy()


# ------------------------------------------------------------
# Construct matching keys
# KEY5: normal matching key
# KEY7: strict key for resolving duplicated KEY5 records
# ------------------------------------------------------------

def make_key5(x):
    return (
        x["CHROM"].astype(str) + ":" +
        x["POS"].astype(str) + ":" +
        x["END"].astype(str) + ":" +
        x["SVTYPE"].astype(str) + ":" +
        x["SVLEN"].astype(str)
    )


def make_key7(x):
    return (
        make_key5(x) + ":" +
        x["REF"].astype(str) + ":" +
        x["ALT"].astype(str)
    )


shared["KEY5"] = make_key5(shared)
shared["KEY7"] = make_key7(shared)


print("Total pan-SVs:", len(df))
print("India-China assembly-shared SVs:", len(shared))
print("Unique KEY5:", shared["KEY5"].nunique())
print("Unique KEY7:", shared["KEY7"].nunique())


# Check duplicated KEY5
dup5 = shared[
    shared.duplicated("KEY5", keep=False)
].copy()

dup5.to_csv(
    "assembly_shared_duplicated_KEY5.tsv",
    sep="\t",
    index=False
)

print("Records involved in duplicated KEY5:", len(dup5))


# Full mother table
shared.to_csv(
    "India_China_shared_assembly_full.tsv",
    sep="\t",
    index=False
)


# Summary by Chinese assembly
summary = []

for sample in chinese:
    n = int(shared[sample + "_present"].sum())

    summary.append({
        "Chinese_assembly": sample,
        "Shared_SVs_with_yd10": n,
        "Total_India_China_shared_SVs": len(shared),
        "Proportion": n / len(shared)
    })

pd.DataFrame(summary).to_csv(
    "India_China_shared_by_Chinese_assembly.tsv",
    sep="\t",
    index=False
)


# Breadth summary
breadth = (
    shared.groupby("N_Chinese")
    .size()
    .reset_index(name="N_SVs")
)

breadth["Proportion"] = breadth["N_SVs"] / len(shared)

breadth.to_csv(
    "India_China_shared_breadth.tsv",
    sep="\t",
    index=False
)

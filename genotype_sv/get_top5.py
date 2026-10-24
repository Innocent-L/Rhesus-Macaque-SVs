import pandas as pd

# read
di = pd.read_csv("DI_SV_siteLevel.sorted.tsv", sep="\t")
fst = pd.read_csv("SV_Fst.weir.fst", sep="\t")

# delet
fst = fst[(fst["WEIR_AND_COCKERHAM_FST"].notna()) & (fst["WEIR_AND_COCKERHAM_FST"] >= 0)]

# top5% value
di_cutoff = di["DI_SV"].quantile(0.95)
fst_cutoff = fst["WEIR_AND_COCKERHAM_FST"].quantile(0.95)
print("DI_SV cutoff =", di_cutoff)
print("FST cutoff =", fst_cutoff)

# top5% sites
di_top = di[di["DI_SV"] >= di_cutoff]
fst_top = fst[fst["WEIR_AND_COCKERHAM_FST"] >= fst_cutoff]

# CHROM + POS
merged = pd.merge(di_top, fst_top, on=["CHROM", "POS"], how="inner")
print("Overlap sites:", len(merged))

# input
merged.to_csv("SV_highly_divergent_overlap.tsv", sep="\t", index=False)


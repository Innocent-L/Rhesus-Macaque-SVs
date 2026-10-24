python3 - <<'PY'
import pandas as pd

x = pd.read_csv(
    "BayPass_China46/china46_allpoly.marker.tsv",
    sep="\t"
)

for p, expected in [
    ("cx",24),
    ("fj",20),
    ("hb",14),
    ("hn",16),
    ("zm",18)
]:
    print(
        p,
        "max_AN =",
        x[f"{p}_AN"].max(),
        "expected maximum =",
        expected
    )

print(
    "total max AN =",
    x["TOTAL_AN"].max()
)
PY

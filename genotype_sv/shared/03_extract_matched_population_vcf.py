import pandas as pd
import pysam

table = pd.read_csv(
    "shared_vs_population_QC_matched.tsv",
    sep="\t",
    dtype=str
)

keep = set(
    table["POP_KEY7"]
    .dropna()
    .astype(str)
)


def make_key(rec):

    svtype = rec.info.get("SVTYPE")

    svlen = rec.info.get("SVLEN")

    if isinstance(svlen, tuple):
        svlen = svlen[0]

    # pysam normally stores END as rec.stop
    end = rec.stop

    ref = rec.ref
    alt = ",".join(rec.alts) if rec.alts else "."

    return (
        f"{rec.chrom}:"
        f"{rec.pos}:"
        f"{end}:"
        f"{svtype}:"
        f"{svlen}:"
        f"{ref}:"
        f"{alt}"
    )


infile = "all_samples_rename.AF_filter.vcf"
outfile = "India_China_shared.population_QC.vcf.gz"

vcf = pysam.VariantFile(infile)

out = pysam.VariantFile(
    outfile,
    "wz",
    header=vcf.header
)

n = 0

for rec in vcf:

    key = make_key(rec)

    if key in keep:
        out.write(rec)
        n += 1

out.close()
vcf.close()

print("Extracted population-QC shared SVs:", n)
print("Expected:", len(keep))

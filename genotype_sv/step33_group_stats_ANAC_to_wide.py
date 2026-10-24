#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd

"""
Input (long):
CHROM POS ID SVTYPE GROUP N_TOTAL N_CALLED AN AC AF ...

Output (wide):
CHROM POS ID SVTYPE AC_cx AN_cx AC_hb AN_hb ... (按你指定顺序)
"""

def main():
    if len(sys.argv) < 3:
        print("Usage: python group_stats_long_to_wide.py sv.group_stats.tsv out.wide.tsv [cx,hb,zm,fj,hn,yd]", file=sys.stderr)
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]

    # 你想要的群体顺序（默认）
    if len(sys.argv) >= 4:
        groups = sys.argv[3].split(",")
    else:
        groups = ["cx","hb","zm","fj","hn","yd"]

    df = pd.read_csv(infile, sep="\t", dtype=str)
    # 转数值列
    for col in ["AN", "AC"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        else:
            raise SystemExit(f"ERROR: column '{col}' not found in input.")

    # 关键列检查
    need = ["CHROM","POS","ID","SVTYPE","GROUP","AN","AC"]
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise SystemExit(f"ERROR: missing columns: {miss}")

    # 确保唯一键：一个SV用 CHROM+POS+ID+SVTYPE（如你还要加 END，就再加一列）
    key_cols = ["CHROM","POS","ID","SVTYPE"]

    # 如果同一 SV×GROUP 出现多行（极少数情况下），这里取最大AN/AC（更保守）
    df2 = (df
           .groupby(key_cols + ["GROUP"], as_index=False)[["AN","AC"]]
           .max())

    # 透视：得到 (key) × (GROUP) 的 AC/AN
    ac_w = df2.pivot_table(index=key_cols, columns="GROUP", values="AC", aggfunc="first")
    an_w = df2.pivot_table(index=key_cols, columns="GROUP", values="AN", aggfunc="first")

    # 统一列名：AC_xx / AN_xx
    ac_w.columns = [f"AC_{g}" for g in ac_w.columns]
    an_w.columns = [f"AN_{g}" for g in an_w.columns]

    wide = pd.concat([ac_w, an_w], axis=1).reset_index()

    # 按你想要的顺序排列列
    ordered_cols = key_cols[:]
    for g in groups:
        ordered_cols.append(f"AC_{g}")
        ordered_cols.append(f"AN_{g}")

    # 如果某些群体在数据里不存在，补空列
    for c in ordered_cols:
        if c not in wide.columns:
            wide[c] = 0

    # 把 AC/AN 缺失补 0，并转为整数
    for g in groups:
        for c in (f"AC_{g}", f"AN_{g}"):
            wide[c] = pd.to_numeric(wide[c], errors="coerce").fillna(0).astype(int)

    wide = wide[ordered_cols]

    wide.to_csv(outfile, sep="\t", index=False)
    print(f"Done. Output -> {outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()


#!/usr/bin/env bash
set -euo pipefail

MARKER="china46_allpoly.marker.tsv"
XTX="china46_core_summary_pi_xtx.out"
OUT="china46_core_XtX_with_SV_ID.tsv"

############################################
# 检查文件
############################################

if [[ ! -f "$MARKER" ]]; then
    echo "ERROR: Cannot find $MARKER"
    exit 1
fi

if [[ ! -f "$XTX" ]]; then
    echo "ERROR: Cannot find $XTX"
    exit 1
fi


############################################
# 合并
#
# marker.tsv:
# MRK CHROM POS ID ...
#
# XtX:
# MRK M_P SD_P M_XtX SD_XtX XtXst log10(1/pval)
############################################

awk '
BEGIN {
    OFS="\t"
}

############################################
# 第一份文件：marker metadata
############################################

NR==FNR {

    if (FNR==1) {
        marker_header=$0
        next
    }

    marker[$1]=$0
    n_marker++
    next
}


############################################
# 第二份文件：BayPass XtX
############################################

FNR==1 {

    printf "%s", marker_header

    for(i=2; i<=NF; i++) {
        printf "%s%s", OFS, $i
    }

    printf "\n"

    next
}


{
    mrk=$1

    if (!(mrk in marker)) {
        print "ERROR: MRK " mrk " not found in marker file" > "/dev/stderr"
        exit 1
    }

    printf "%s", marker[mrk]

    for(i=2; i<=NF; i++) {
        printf "%s%s", OFS, $i
    }

    printf "\n"

    n_xtx++
}


END {

    if (n_marker != n_xtx) {
        print "WARNING: marker rows = " n_marker \
              ", XtX rows = " n_xtx > "/dev/stderr"
    }
}
' "$MARKER" "$XTX" > "$OUT"


############################################
# 检查
############################################

echo "Finished."
echo "Output: $OUT"

echo
echo "Marker rows:"
tail -n +2 "$MARKER" | wc -l

echo "XtX rows:"
tail -n +2 "$XTX" | wc -l

echo "Merged rows:"
tail -n +2 "$OUT" | wc -l

echo
echo "Preview:"
head -5 "$OUT"

#!/usr/bin/env bash
set -euo pipefail

g_baypass \
    -gfile 04_pod/G.china46_POD \
    -outprefix 04_pod/china46_POD \
    -seed 20260929 \
    -nthreads 40 \
    -npilot 20 \
    -pilotlength 500 \
    -burnin 5000 \
    -thin 25 \
    -nval 1000 \
    > 04_pod/china46_POD.log 2>&1

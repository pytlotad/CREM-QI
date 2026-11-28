#!/bin/bash
# Audit 873 POST HOC (after the prereg verdicts): X arm was censored in orbit count (n = 2 period 8x longer); rerun X with TEND x5 (~5800 orbits).
set -u; export LC_ALL=C; ulimit -v 4000000
D=$(cd "$(dirname "$0")/.." && pwd); B=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s873/local871; O=$D/out/runs; export B O
PS="M1=1 G1=2.00231930436"
{ echo "XL_H_2.02 N0=2.02 E0=0.95 W=0.5 DX=0.1 TEND=300000"; echo "XL_P_2.02 $PS N0=2.02 E0=0.95 W=1.0 DX=0.2 TEND=600000"; } |
xargs -P 2 -L 1 bash -c 'ulimit -v 2500000; name=$0; env "$@" $B > $O/$name.txt 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

#!/bin/bash
# Audit 833: static comparison of the two Thomas expressions (p833 from the tree; engine unchanged).
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s833; D=$R/audits/data/section-833; O=$D/out/runs
mkdir -p $S $O
( cd $R && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-833/p833.cpp -o $S/p833 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"; export S O
C=0.769637161900118
N5="u5a,0.5120091953982716 u5b,0.5652161499007791 u5c,0.66666666666666663 u5d,0.8125025652578265 u5e,0.9551918755393795"
{ for nl in $N5; do IFS=, read nn f <<< "$nl"
    echo "ops_$nn electron,positron -$C -$C $f"; echo "pps_$nn electron,positron -$C $C $f"
    echo "H0_$nn proton,electron $C -$C $f"; echo "H1_$nn proton,electron $C $C $f"; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; $S/p833 "$@" $name >$O/$name.out 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

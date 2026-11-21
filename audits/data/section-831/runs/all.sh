#!/bin/bash
# Audit 831: J4L3 payout weight with a slower window (eps 0.002) on the mechanical path; p787 from the tree.  Waits for 830 ALLDONE.
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s831; D=$R/audits/data/section-831; O=$D/out/runs
until grep -q ALLDONE $R/audits/data/section-830/out/runs830.log; do sleep 60; done
echo "830 done $(date +%T)"
rm -rf $S; mkdir -p $S $O/WS
( cd $R && g++ -O3 -march=native -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-787/p787.cpp -o $S/p787 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"; export S O
{ for sys in "H0 PAIR=proton,electron PHEN=1" "H1 PAIR=proton,electron PHEN=2" "pPs PHEN=1" "oPs PHEN=2"; do set -- $sys; tag=$1; shift
    for i in 0 1 2; do echo "WS/${tag}_4f_$i $i CREM_WINDOW_WEIGHT=J4L3 CREM_LADDER_WINDOW_FRACTION=0.002 $*"; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env CREM_MECHANICAL_HAZARD_SCALE=8e5 CREM_DEBUG_LADDER=1 CREM_DEBUG_LADDER_SPLIT=1 CREM_PHASE_TRACE=1 NLEV=4 LHBAR=3 OBS=30 "$@" $S/p787 $seed >$O/$name.out 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

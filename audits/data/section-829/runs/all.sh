#!/bin/bash
# Audit 829: window payout weight (CREM_WINDOW_WEIGHT=J4L3 / phi) on the mechanical path; p787 from the tree.  Waits for 828 ALLDONE.
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s829; D=$R/audits/data/section-829; O=$D/out/runs
until grep -q ALLDONE $R/audits/data/section-828/out/runs828.log; do sleep 60; done
echo "828 done $(date +%T)"
rm -rf $S; mkdir -p $S $O/WJL $O/WF $O/W0
( cd $R && g++ -O3 -march=native -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-787/p787.cpp -o $S/p787 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"; export S O
{ echo "W0/H1_4f_0 0 PAIR=proton,electron PHEN=2"; echo "W0/oPs_4f_0 0 PHEN=2"
  for v in "WJL J4L3" "WF phi"; do set -- $v; V=$1; X=CREM_WINDOW_WEIGHT=$2
  for sys in "H0 PAIR=proton,electron PHEN=1" "H1 PAIR=proton,electron PHEN=2" "pPs PHEN=1" "oPs PHEN=2"; do set -- $sys; tag=$1; shift
    for i in 0 1 2; do echo "$V/${tag}_4f_$i $i $X $*"; done; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env CREM_MECHANICAL_HAZARD_SCALE=8e5 CREM_DEBUG_LADDER=1 CREM_DEBUG_LADDER_SPLIT=1 CREM_PHASE_TRACE=1 NLEV=4 LHBAR=3 OBS=30 "$@" $S/p787 $seed >$O/$name.out 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

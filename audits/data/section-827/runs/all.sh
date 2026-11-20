#!/bin/bash
# Audit 827: source of the hydrogen drift in the payout window.  p787 from the tree (code 826); V1 tol 1e-7, V2 512 steps/orbit,
# V3 window fraction 0.002; V0 = 4f runs of 826.  36 runs, xargs -P 3.
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s827; D=$R/audits/data/section-827; O=$D/out/runs
rm -rf $S; mkdir -p $S $O/V1 $O/V2 $O/V3
( cd $R && g++ -O3 -march=native -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-787/p787.cpp -o $S/p787 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"; export S O
{ for v in "V1 CREM_MECHANICAL_TOLERANCE=1e-7" "V2 CREM_MECHANICAL_STEPS_PER_ORBIT=512" "V3 CREM_LADDER_WINDOW_FRACTION=0.002"; do set -- $v; V=$1; X=$2
  for sys in "H0 PAIR=proton,electron PHEN=1" "H1 PAIR=proton,electron PHEN=2" "pPs PHEN=1" "oPs PHEN=2"; do set -- $sys; tag=$1; shift
    for i in 0 1 2; do echo "$V/${tag}_4f_$i $i NLEV=4 LHBAR=3 OBS=30 $X $*"; done; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env CREM_MECHANICAL_HAZARD_SCALE=8e5 CREM_DEBUG_LADDER=1 CREM_DEBUG_LADDER_SPLIT=1 CREM_PHASE_TRACE=1 "$@" $S/p787 $seed >$O/$name.out 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

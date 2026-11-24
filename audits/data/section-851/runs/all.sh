#!/bin/bash
# Audit 851: CREM_PHOTON_INCREMENT_STABLE (dynamics switch, default off).  S: static probe p851 (H, muonium, Ps; with / without).
# W: as 850 (eps 0.001, start 4f, DRIFT850) + switch, H F0/F1, muonium F0/F1, p-Ps, o-Ps, seeds 0-9; K = 850 W/H0_4f_0 without it.
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s851; D=$R/audits/data/section-851; O=$D/out/runs
mkdir -p $S $O/K $O/W $O/S
( cd $R && g++ -O3 -march=native -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-787/p787.cpp -o $S/p787 $(root-config --libs) \
  && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules $D/p851.cpp -o $S/p851 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"; export S O
for p in proton,electron antimuon,electron electron,positron; do
  ( ulimit -v 2500000; $S/p851 $p 10000 ) > $O/S/${p/,/_}_0.txt 2>&1
  ( ulimit -v 2500000; CREM_PHOTON_INCREMENT_STABLE=1 $S/p851 $p 10000 ) > $O/S/${p/,/_}_1.txt 2>&1
done; echo "static $(date +%T)"
SYS=("H0 PAIR=proton,electron PHEN=1" "H1 PAIR=proton,electron PHEN=2" "M0 PAIR=antimuon,electron PHEN=1" "M1 PAIR=antimuon,electron PHEN=2" "pPs PHEN=1" "oPs PHEN=2")
{ echo "K/H0_4f_0 0 NLEV=4 LHBAR=3 OBS=30 CREM_LADDER_WINDOW_FRACTION=0.001 CREM_DEBUG_WINDOW_DRIFT=1 PAIR=proton,electron PHEN=1"
  for sys in "${SYS[@]}"; do set -- $sys; tag=$1; shift
    for i in 0 1 2 3 4 5 6 7 8 9; do echo "W/${tag}_4f_$i $i NLEV=4 LHBAR=3 OBS=30 CREM_LADDER_WINDOW_FRACTION=0.001 CREM_DEBUG_WINDOW_DRIFT=1 CREM_PHOTON_INCREMENT_STABLE=1 $*"; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env CREM_MECHANICAL_HAZARD_SCALE=8e5 CREM_DEBUG_LADDER=1 CREM_DEBUG_LADDER_SPLIT=1 CREM_PHASE_TRACE=1 "$@" $S/p787 $seed >$O/$name.out 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

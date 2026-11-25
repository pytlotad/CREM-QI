#!/bin/bash
# Audit 852 M: mechanical path with the new default (H0, H1, M0, oPs seed 0) vs 851; L = H0 with CREM_PHOTON_INCREMENT_LEGACY=1 vs 850.
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s852; D=$R/audits/data/section-852; O=$D/out/runs
mkdir -p $S $O/M
( cd $R && g++ -O3 -march=native -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-787/p787.cpp -o $S/p787 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"; export S O
{ echo "M/H0_4f_0 0 PAIR=proton,electron PHEN=1"; echo "M/H1_4f_0 0 PAIR=proton,electron PHEN=2"; echo "M/M0_4f_0 0 PAIR=antimuon,electron PHEN=1"
  echo "M/oPs_4f_0 0 PHEN=2"; echo "M/L_H0_4f_0 0 PAIR=proton,electron PHEN=1 CREM_PHOTON_INCREMENT_LEGACY=1"; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env NLEV=4 LHBAR=3 OBS=30 CREM_LADDER_WINDOW_FRACTION=0.001 CREM_DEBUG_WINDOW_DRIFT=1 CREM_MECHANICAL_HAZARD_SCALE=8e5 CREM_DEBUG_LADDER=1 CREM_DEBUG_LADDER_SPLIT=1 CREM_PHASE_TRACE=1 "$@" $S/p787 $seed >$O/$name.out 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "MDONE $(date +%T)"

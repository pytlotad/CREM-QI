#!/bin/bash
# Audit 849: payout-window drift vs window fraction eps (0.01 / 0.005 / 0.002 / 0.001), start 4f, H F0/F1, p-Ps, o-Ps, muonium F0/F1.
# H / Ps at 0.01 and 0.002 from the archives 826 / 827 (K1, K2 = byte checks).  p787 from the tree (as 827).
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s849; D=$R/audits/data/section-849; O=$D/out/runs
rm -rf $S; mkdir -p $S $O/K $O/e0.01 $O/e0.005 $O/e0.002 $O/e0.001
( cd $R && g++ -O3 -march=native -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-787/p787.cpp -o $S/p787 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"; export S O
SYS=("H0 PAIR=proton,electron PHEN=1" "H1 PAIR=proton,electron PHEN=2" "pPs PHEN=1" "oPs PHEN=2" "M0 PAIR=antimuon,electron PHEN=1" "M1 PAIR=antimuon,electron PHEN=2")
{ echo "K/H0_4f_0_e0.01 0 NLEV=4 LHBAR=3 OBS=30 PAIR=proton,electron PHEN=1"
  echo "K/H0_4f_0_e0.002 0 NLEV=4 LHBAR=3 OBS=30 CREM_LADDER_WINDOW_FRACTION=0.002 PAIR=proton,electron PHEN=1"
  for e in 0.001 0.005 0.002 0.01; do
    for sys in "${SYS[@]}"; do set -- $sys; tag=$1; shift
      case $tag in M*) ;; *) [ $e = 0.01 -o $e = 0.002 ] && continue;; esac
      X=""; [ $e != 0.01 ] && X=CREM_LADDER_WINDOW_FRACTION=$e
      for i in 0 1 2; do echo "e$e/${tag}_4f_$i $i NLEV=4 LHBAR=3 OBS=30 $X $*"; done; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env CREM_MECHANICAL_HAZARD_SCALE=8e5 CREM_DEBUG_LADDER=1 CREM_DEBUG_LADDER_SPLIT=1 CREM_PHASE_TRACE=1 "$@" $S/p787 $seed >$O/$name.out 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

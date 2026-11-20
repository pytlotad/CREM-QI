#!/bin/bash
# Audit 825: mechanical-path cascade landing diagnosis (print only).  p787 from the tree; 48 runs, xargs -P 3.
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s825; D=$R/audits/data/section-825; O=$D/out/runs
mkdir -p $S $O
( cd $R && g++ -O3 -march=native -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-787/p787.cpp -o $S/p787 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"; export S O
{ for sys in "H0 PAIR=proton,electron PHEN=1" "H1 PAIR=proton,electron PHEN=2" "pPs PHEN=1" "oPs PHEN=2"; do set -- $sys; tag=$1; shift
    for st in "4f 4 3 30" "4d 4 2 30" "4p 4 1 30" "3d 3 2 40"; do set -- $st $*; s=$1 n=$2 l=$3 obs=$4; shift 4
      for i in 0 1 2; do echo "${tag}_${s}_$i $i NLEV=$n LHBAR=$l OBS=$obs $*"; done; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env CREM_MECHANICAL_HAZARD_SCALE=8e5 CREM_DEBUG_LADDER=1 CREM_DEBUG_LADDER_SPLIT=1 CREM_PHASE_TRACE=1 "$@" $S/p787 $seed >$O/$name.out 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

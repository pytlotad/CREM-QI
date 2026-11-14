#!/bin/bash
# Audit 787 runs: probe p787 (= ../section-624/p624.cpp + PHEN) from the working tree (p787) and from HEAD 9f8537e6 (p787_head, control F).
# 3d start (NLEV 3, LHBAR 2), OBS 40 orbits n = 3, hazard scale 8e5, trace CREM_PHASE_TRACE.  P: + CREM_PHASE_THRESHOLD (200 seeds per system), D: trace only (60),
# F: no switches, p787_head vs p787 (H, seeds 0-9), R: P repeated (H, seeds 0-9).  Systems: H (PAIR=proton,electron), pPs (PHEN 1), oPs (PHEN 2).  3 parallel, 2.5 GB, LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r787; O=$D/out/runs; mkdir -p $O/P $O/D $O/F $O/R
export S O
echo "start $(date +%T)"
{ for sys in "H PAIR=proton,electron" "pPs PHEN=1" "oPs PHEN=2"; do set -- $sys
    for i in $(seq 0 199); do echo "P/${1}_$i p787 $i $2 CREM_PHASE_TRACE=1 CREM_PHASE_THRESHOLD=1"; done
    for i in $(seq 0 59); do echo "D/${1}_$i p787 $i $2 CREM_PHASE_TRACE=1"; done; done
  for i in $(seq 0 9); do echo "F/head_$i p787_head $i PAIR=proton,electron"; echo "F/new_$i p787 $i PAIR=proton,electron"
    echo "R/H_$i p787 $i PAIR=proton,electron CREM_PHASE_TRACE=1 CREM_PHASE_THRESHOLD=1"; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; bin=$1; seed=$2; shift 2
  env NLEV=3 LHBAR=2 OBS=40 CREM_MECHANICAL_HAZARD_SCALE=8e5 "$@" $S/$bin $seed >$O/$name.out 2>$O/$name.err && echo "done $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

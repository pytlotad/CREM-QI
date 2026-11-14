#!/bin/bash
# Audit 787, reduced set (amendment 2): P 100 seeds per system, D 30, F 10 x 2, R 10; runs already finished (out/runs787.log "done") are skipped.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r787; O=$D/out/runs; mkdir -p $O/P $O/D $O/F $O/R
export S O; DONE=$(grep '^done ' $D/out/runs787.log | awk '{print $2}')
echo "start2 $(date +%T)"
{ for sys in "H PAIR=proton,electron" "pPs PHEN=1" "oPs PHEN=2"; do set -- $sys
    for i in $(seq 0 99); do echo "P/${1}_$i p787 $i $2 CREM_PHASE_TRACE=1 CREM_PHASE_THRESHOLD=1"; done
    for i in $(seq 0 29); do echo "D/${1}_$i p787 $i $2 CREM_PHASE_TRACE=1"; done; done
  for i in $(seq 0 9); do echo "F/head_$i p787_head $i PAIR=proton,electron"; echo "F/new_$i p787 $i PAIR=proton,electron"
    echo "R/H_$i p787 $i PAIR=proton,electron CREM_PHASE_TRACE=1 CREM_PHASE_THRESHOLD=1"; done; } |
awk -v done="$DONE" 'BEGIN{n=split(done,a,"\n"); for(i=1;i<=n;i++) d[a[i]]=1} !($1 in d)' |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; bin=$1; seed=$2; shift 2
  env NLEV=3 LHBAR=2 OBS=40 CREM_MECHANICAL_HAZARD_SCALE=8e5 "$@" $S/$bin $seed >$O/$name.out 2>$O/$name.err && echo "done $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

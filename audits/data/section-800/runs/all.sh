#!/bin/bash
# Audit 800: rule on the per-orbit measure (799) with more seeds: B (F = 1, seeds 7-12), B0 (F = 0, seeds 4 and 6).  Probe p799 (tree ac3777ab), no code change.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r799; O=$D/out/runs
mkdir -p $O/B $O/B0; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_MECHANICAL_STEPS_PER_ORBIT=512"
TORI="3,L1,0.3333333333333333 3,L1.5,0.5 3,L2,0.6666666666666666 3,L2.5,0.8333333333333334 3,L3,1 2,L1,0.5 2,L1.5,0.75 2,L2,1"
TYP=("proton PAIR=proton,electron" "antiproton PAIR=antiproton,positron")
echo "start $(date +%T)"
{ for s in 7 8 9 10 11 12; do for nl in $TORI; do IFS=, read n L f <<< "$nl"
    for t in "${TYP[@]}"; do set -- $t; echo "B/${1}_n$n${L}_$s $s $B NLEV=$n LFRAC=$f PHEN=2 $2"; done; done; done
  for s in 4 6; do for nl in $TORI; do IFS=, read n L f <<< "$nl"
    for t in "${TYP[@]}"; do set -- $t; echo "B0/${1}_n$n${L}_$s $s $B NLEV=$n LFRAC=$f PHEN=1 $2"; done; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p799 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

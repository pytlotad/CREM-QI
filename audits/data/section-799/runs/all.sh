#!/bin/bash
# Audit 799: precession per radial orbit (CREM_PAIR_TORQUE_PASS=1, with CREM_PAIR_TORQUE_APO=0.5): A (F = 1, seeds 4-6), A0 (F = 0, seed 5),
# AZ (n3L1 seed 5 at 2048 steps), AK (old start CREM_KEPLER_START=1, n3L1 and n2L1, seeds 4-6).  Probe p799 from the working tree (audit 799 switch); same runs as 798 without K.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r799; O=$D/out/runs
mkdir -p $O/A $O/A0 $O/AZ $O/AK; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1"
TORI="3,L1,0.3333333333333333 3,L1.5,0.5 3,L2,0.6666666666666666 3,L2.5,0.8333333333333334 3,L3,1 2,L1,0.5 2,L1.5,0.75 2,L2,1"
TYP=("proton PAIR=proton,electron" "antiproton PAIR=antiproton,positron")
echo "start $(date +%T)"
{ for t in "${TYP[@]}"; do set -- $t
    echo "AZ/${1}_n3L1_5 5 $B CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_MECHANICAL_STEPS_PER_ORBIT=2048 NLEV=3 LFRAC=0.3333333333333333 PHEN=2 $2"; done
  for s in 4 5 6; do for nl in $TORI; do IFS=, read n L f <<< "$nl"
    for t in "${TYP[@]}"; do set -- $t
      echo "A/${1}_n$n${L}_$s $s $B CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_MECHANICAL_STEPS_PER_ORBIT=512 NLEV=$n LFRAC=$f PHEN=2 $2"; done; done; done
  for s in 4 5 6; do for nl in "3,L1,0.3333333333333333" "2,L1,0.5"; do IFS=, read n L f <<< "$nl"
    for t in "${TYP[@]}"; do set -- $t
      echo "AK/${1}_n$n${L}_$s $s $B CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_KEPLER_START=1 CREM_MECHANICAL_STEPS_PER_ORBIT=512 NLEV=$n LFRAC=$f PHEN=2 $2"; done; done; done
  for nl in $TORI; do IFS=, read n L f <<< "$nl"
    for t in "${TYP[@]}"; do set -- $t
      echo "A0/${1}_n$n${L}_5 5 $B CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_MECHANICAL_STEPS_PER_ORBIT=512 NLEV=$n LFRAC=$f PHEN=1 $2"; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p799 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

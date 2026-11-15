#!/bin/bash
# Audit 794: A (old start, F = 1, hydrogen seed 5, trace), B (CREM_PAIR_TORQUE_RELSTART, F = 1), C (RELSTART, F = 0); all CREM_NO_DIPOLE_FORCE=1, 512 steps.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r794; O=$D/out/runs
mkdir -p $O/A $O/B $O/C; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_NO_DIPOLE_FORCE=1 CREM_MECHANICAL_STEPS_PER_ORBIT=512"
TR="CREM_PAIR_TORQUE_CHUNK=1 CREM_PAIR_TORQUE_ENERGY=1"; RS="CREM_PAIR_TORQUE_RELSTART=1"
TORI="3,L1,0.3333333333333333 3,L1.5,0.5 3,L2,0.6666666666666666 3,L2.5,0.8333333333333334 3,L3,1 2,L1,0.5 2,L1.5,0.75 2,L2,1"
echo "start $(date +%T)"
{ for nl in $TORI; do IFS=, read n L f <<< "$nl"
    echo "A/proton_n$n${L}_5 5 $B NLEV=$n LFRAC=$f PHEN=2 PAIR=proton,electron $TR"
    echo "B/proton_n$n${L}_5 5 $B NLEV=$n LFRAC=$f PHEN=2 PAIR=proton,electron $RS $TR"; done
  for nl in $TORI; do IFS=, read n L f <<< "$nl"
    for s in 4 5 6; do for t in "proton PAIR=proton,electron" "antiproton PAIR=antiproton,positron"; do set -- $t
      [ "$1_$s" != "proton_5" ] && echo "B/${1}_n$n${L}_$s $s $B NLEV=$n LFRAC=$f PHEN=2 $2 $RS"
      echo "C/${1}_n$n${L}_$s $s $B NLEV=$n LFRAC=$f PHEN=1 $2 $RS"; done; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p794 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

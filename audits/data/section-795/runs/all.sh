#!/bin/bash
# Audit 795: S (start line TORQUE794, hydrogen F = 1 / F = 0, production start), S0 (CREM_KEPLER_START=1, F = 1), G (rule, production start, F = 1, seeds 4-6),
# G0 (F = 0, seed 5).  Unit check t795 -> out/t795.txt.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r795; O=$D/out/runs
mkdir -p $O/S $O/S0 $O/G $O/G0; export S O
$S/t795 > $D/out/t795.txt
C="MODE=est QSPIN=1 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1"
SB="$C BUDGET=5 CREM_PAIR_TORQUE=1 CREM_PAIR_TORQUE_RELSTART=1"; GB="$C BUDGET=60 CREM_PAIR_TORQUE=200 CREM_NO_DIPOLE_FORCE=1 CREM_MECHANICAL_STEPS_PER_ORBIT=512"
TORI="3,L1,0.3333333333333333 3,L1.5,0.5 3,L2,0.6666666666666666 3,L2.5,0.8333333333333334 3,L3,1 2,L1,0.5 2,L1.5,0.75 2,L2,1"
echo "start $(date +%T)"
{ for nl in $TORI; do IFS=, read n L f <<< "$nl"
    echo "S/F1_n$n${L} 5 $SB NLEV=$n LFRAC=$f PHEN=2 PAIR=proton,electron"
    echo "S/F0_n$n${L} 5 $SB NLEV=$n LFRAC=$f PHEN=1 PAIR=proton,electron"
    echo "S0/F1_n$n${L} 5 $SB NLEV=$n LFRAC=$f PHEN=2 PAIR=proton,electron CREM_KEPLER_START=1"; done
  for nl in $TORI; do IFS=, read n L f <<< "$nl"
    for s in 4 5 6; do for t in "proton PAIR=proton,electron" "antiproton PAIR=antiproton,positron"; do set -- $t
      echo "G/${1}_n$n${L}_$s $s $GB NLEV=$n LFRAC=$f PHEN=2 $2"
      [ $s = 5 ] && echo "G0/${1}_n$n${L}_$s $s $GB NLEV=$n LFRAC=$f PHEN=1 $2"; done; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p795 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

#!/bin/bash
# Audit 796: unit check t796 -> out/t796.txt; G (rule, production start, hydrogen / antihydrogen, F = 1, seed 5, 8 tori) -> out/runs/G.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r796; O=$D/out/runs
mkdir -p $O/G; export S O
$S/t796 > $D/out/t796.txt
C="MODE=est QSPIN=1 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1"
GB="$C BUDGET=60 CREM_PAIR_TORQUE=200 CREM_NO_DIPOLE_FORCE=1 CREM_MECHANICAL_STEPS_PER_ORBIT=512"
TORI="3,L1,0.3333333333333333 3,L1.5,0.5 3,L2,0.6666666666666666 3,L2.5,0.8333333333333334 3,L3,1 2,L1,0.5 2,L1.5,0.75 2,L2,1"
echo "start $(date +%T)"
{ for nl in $TORI; do IFS=, read n L f <<< "$nl"
    for t in "proton PAIR=proton,electron" "antiproton PAIR=antiproton,positron"; do set -- $t
      echo "G/${1}_n$n${L}_5 5 $GB NLEV=$n LFRAC=$f PHEN=2 $2"; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p796 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

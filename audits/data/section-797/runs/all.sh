#!/bin/bash
# Audit 797: rule WITH moment forces on the new start (795 + 796): K (CREM_KEPLER_START=1 control vs 791), G (F = 1, seeds 4-6), G0 (F = 0, seed 5),
# Z (n3L1 seed 5 at 2048 steps).  Probe p796 (unchanged) from tree f8b0ccf9.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r796; O=$D/out/runs
mkdir -p $O/K $O/G $O/G0 $O/Z; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1"
TORI="3,L1,0.3333333333333333 3,L1.5,0.5 3,L2,0.6666666666666666 3,L2.5,0.8333333333333334 3,L3,1 2,L1,0.5 2,L1.5,0.75 2,L2,1"
echo "start $(date +%T)"
{ echo "K/proton_n3L1_5 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=512 NLEV=3 LFRAC=0.3333333333333333 PHEN=2 PAIR=proton,electron CREM_KEPLER_START=1"
  for t in "proton PAIR=proton,electron" "antiproton PAIR=antiproton,positron"; do set -- $t
    echo "Z/${1}_n3L1_5 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=2048 NLEV=3 LFRAC=0.3333333333333333 PHEN=2 $2"; done
  for s in 4 5 6; do for nl in $TORI; do IFS=, read n L f <<< "$nl"
    for t in "proton PAIR=proton,electron" "antiproton PAIR=antiproton,positron"; do set -- $t
      echo "G/${1}_n$n${L}_$s $s $B CREM_MECHANICAL_STEPS_PER_ORBIT=512 NLEV=$n LFRAC=$f PHEN=2 $2"; done; done; done
  for nl in $TORI; do IFS=, read n L f <<< "$nl"
    for t in "proton PAIR=proton,electron" "antiproton PAIR=antiproton,positron"; do set -- $t
      echo "G0/${1}_n$n${L}_5 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=512 NLEV=$n LFRAC=$f PHEN=1 $2"; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p796 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

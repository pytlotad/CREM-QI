#!/bin/bash
# Audit 791: F (2 file-control runs), Z (3 convergence runs at 2048 steps), G (48 main runs at 512 steps) -- CREM_PAIR_TORQUE mechB0 lines.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r791; O=$D/out/runs
mkdir -p $O/F $O/Z $O/G; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 PHEN=2"
echo "start $(date +%T)"
{ echo "F/proton_n2L1.5_2_off 2 $B NLEV=2 LFRAC=0.75 PAIR=proton,electron"
  echo "F/proton_n2L1.5_2_on 2 $B NLEV=2 LFRAC=0.75 PAIR=proton,electron CREM_PAIR_TORQUE_ANYL=1"
  for nl in "L2 0.6666666666666666" "L2.5 0.8333333333333334" "L3 1"; do set -- $nl
    echo "Z/proton_n3$1_4 4 $B NLEV=3 LFRAC=$2 PAIR=proton,electron CREM_PAIR_TORQUE_ANYL=1 CREM_MECHANICAL_STEPS_PER_ORBIT=2048"; done
  for nl in "3 L1 0.3333333333333333" "3 L1.5 0.5" "3 L2 0.6666666666666666" "3 L2.5 0.8333333333333334" "3 L3 1" "2 L1 0.5" "2 L1.5 0.75" "2 L2 1"; do
    for s in 4 5 6; do for t in "proton PAIR=proton,electron" "antiproton PAIR=antiproton,positron"; do set -- $t $nl
      echo "G/${1}_n$3$4_$s $s $B NLEV=$3 LFRAC=$5 $2 CREM_PAIR_TORQUE_ANYL=1 CREM_MECHANICAL_STEPS_PER_ORBIT=512"; done; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p791 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

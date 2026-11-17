#!/bin/bash
# Audit 802: secular-path raw averages (Kepler / relativistic / Sommerfeld nodes, CREM_SECULAR_NODES=1) at the 699 start.  K: control (N = 200); S: short (N = 2).
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r802; O=$D/out/runs
mkdir -p $O/K $O/S; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_MECHANICAL_STEPS_PER_ORBIT=512 CREM_SECULAR_NODES=1"
TORI="3,L1,0.3333333333333333 3,L1.5,0.5 3,L2,0.6666666666666666 3,L2.5,0.8333333333333334 3,L3,1 2,L1,0.5 2,L1.5,0.75 2,L2,1"
TYP=("proton PAIR=proton,electron" "antiproton PAIR=antiproton,positron")
echo "start $(date +%T)"
{ echo "K/proton_n3L1_5 5 $B CREM_PAIR_TORQUE=200 NLEV=3 LFRAC=0.3333333333333333 PHEN=2 PAIR=proton,electron"
  for F in 1 0; do if [ $F = 1 ]; then SEEDS="4 5 6 7 8 9 10 11 12"; PH=2; else SEEDS="4 5 6"; PH=1; fi
  for s in $SEEDS; do for nl in $TORI; do IFS=, read n L f <<< "$nl"
    for t in "${TYP[@]}"; do set -- $t; echo "S/F${F}_${1}_n$n${L}_$s $s $B CREM_PAIR_TORQUE=2 NLEV=$n LFRAC=$f PHEN=$PH $2"; done; done; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p802 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

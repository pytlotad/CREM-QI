#!/bin/bash
# Audit 807: rule D vs Born B for hydrogen 2P at 2048 steps, n = 2 tori, more seeds; probe p804.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r804; O=$D/out/runs
mkdir -p $O/V3; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_SECULAR_NODES=1"
TORI="2,L1,0.5 2,L1.5,0.75 2,L2,1"
TYP=("proton PAIR=proton,electron" "antiproton PAIR=antiproton,positron")
echo "start $(date +%T)"
{ for F in 1 0; do PH=$((F + 1)); if [ $F = 1 ]; then SEEDS="4 6 7 8"; else SEEDS="4 6"; fi
  for s in $SEEDS; do for nl in $TORI; do IFS=, read n L f <<< "$nl"
    for t in "${TYP[@]}"; do set -- $t; echo "V3/F${F}_${1}_n$n${L}_$s $s $B CREM_MECHANICAL_STEPS_PER_ORBIT=2048 NLEV=$n LFRAC=$f PHEN=$PH $2"; done; done; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p804 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

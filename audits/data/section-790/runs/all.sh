#!/bin/bash
# Audit 790: 32 torque runs (CREM_PAIR_TORQUE=200, mechB0) -- hydrogen / antihydrogen, seeds 2 and 3, tori n = 2 (L = 1, 1.5, 2) and n = 3 (L = 1 ... 3).
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r789; O=$D/out/runs
mkdir -p $O; export S O
HY="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 PHEN=2"
echo "start $(date +%T)"
{ for s in 2 3; do for t in "proton PAIR=proton,electron" "antiproton PAIR=antiproton,positron"; do set -- $t
    for nl in "3 L1 0.3333333333333333" "3 L1.5 0.5" "3 L2 0.6666666666666666" "3 L2.5 0.8333333333333334" "3 L3 1" "2 L1 0.5" "2 L1.5 0.75" "2 L2 1"; do
      set -- $t $nl; echo "${1}_n$3$4_$s $s $HY NLEV=$3 LFRAC=$5 $2"; done; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p789 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

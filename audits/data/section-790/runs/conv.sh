#!/bin/bash
# Audit 790, POST HOC: step-density convergence of Omega(n, L) -- hydrogen seed 2, CREM_MECHANICAL_STEPS_PER_ORBIT 512 / 2048 (default 128).
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r789; O=$D/out/conv
mkdir -p $O; export S O
HY="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 PHEN=2 PAIR=proton,electron"
echo "start $(date +%T)"
{ for st in 512 2048; do for nl in "2 L1 0.5" "2 L1.5 0.75" "2 L2 1" "3 L1 0.3333333333333333" "3 L1.5 0.5"; do set -- $nl
    echo "proton_n$1$2_2_s$st 2 $HY NLEV=$1 LFRAC=$3 CREM_MECHANICAL_STEPS_PER_ORBIT=$st"; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p789 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

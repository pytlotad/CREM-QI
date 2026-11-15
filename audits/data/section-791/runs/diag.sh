#!/bin/bash
# Audit 791 POST HOC diagnostic (after the scored results): n = 3, L = 1, seed 5, proton and antiproton --
# window 400 orbits at 512 steps, and window 200 orbits at 2048 steps.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r791; O=$D/out/runs
mkdir -p $O/D; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 PHEN=2 NLEV=3 LFRAC=0.3333333333333333 CREM_PAIR_TORQUE_ANYL=1"
echo "start $(date +%T)"
{ for t in "proton PAIR=proton,electron" "antiproton PAIR=antiproton,positron"; do set -- $t
    echo "D/${1}_n3L1_5_T400 5 $B $2 CREM_PAIR_TORQUE=400 CREM_MECHANICAL_STEPS_PER_ORBIT=512"
    echo "D/${1}_n3L1_5_S2048 5 $B $2 CREM_PAIR_TORQUE=200 CREM_MECHANICAL_STEPS_PER_ORBIT=2048"; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p791 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

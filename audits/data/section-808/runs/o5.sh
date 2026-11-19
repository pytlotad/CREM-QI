#!/bin/bash
# Audit 808, POPRAWKA (pilot: integrator order 2): O5 = the 3 Gauss nodes at 4096 steps, F = 1 hydrogen, seed 5; probe p808.  Runs beside all.sh (-P 1).
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r808; O=$D/out/runs
mkdir -p $O/O5; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_SECULAR_NODES=1 NLEV=2 CREM_NO_DIPOLE_FORCE=1 CREM_PAIR_TORQUE_APSE=1"
echo "start $(date +%T)"
for nl in u3a,0.52985793589488495 u3b,0.66666666666666663 u3c,0.89871349267654366; do IFS=, read nn f <<< "$nl"
  echo "O5/F1_proton_${nn}_s4096_5 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=4096 LFRAC=$f PHEN=2 PAIR=proton,electron"; done |
xargs -P 1 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p808 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

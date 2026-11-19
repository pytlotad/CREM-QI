#!/bin/bash
# Audit 809: quadrature control at the step limit (4 Gauss nodes in u = 1/L at 1024 / 4096, F = 1; 2048 F = 0) and the C image of the
# antihydrogen anomaly runs (hydrogen with the heavy spin reversed, CREM_PAIR_TORQUE_HEAVY_FLIP=1); probe p809.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r809; O=$D/out/runs
mkdir -p $O/Q1 $O/Q2 $O/C1; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_SECULAR_NODES=1 NLEV=2 CREM_NO_DIPOLE_FORCE=1 CREM_PAIR_TORQUE_APSE=1"
N4="u4a,0.51798223077348682 u4b,0.59880579377581356 u4c,0.75187434103678785 u4d,0.93507595217092121"
job() { echo "$1/F$2_proton_$3_s$4_5 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=$4 LFRAC=$5 PHEN=$(($2 + 1)) PAIR=proton,electron${6:+ $6}"; }
echo "start $(date +%T)"
{ for nl in $N4; do IFS=, read nn f <<< "$nl"; job Q1 1 $nn 4096 $f; done
  for nl in $N4; do IFS=, read nn f <<< "$nl"; job Q2 0 $nn 2048 $f; done
  for F in 1 0; do job C1 $F flipgE 2048 0.75 "CREM_PAIR_TORQUE_HEAVY_FLIP=1"; job C1 $F flipg2 2048 0.75 "CREM_PAIR_TORQUE_HEAVY_FLIP=1 CREM_PAIR_TORQUE_BMT_G=2"; done
  for nl in $N4; do IFS=, read nn f <<< "$nl"; job Q1 1 $nn 1024 $f; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p809 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

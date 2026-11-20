#!/bin/bash
# Audit 815: Ps orbit under the production force with the Plummer floor scaled (CREM_SEPARATION_FLOOR_SCALE); probe p814 (tree 72b7c76c).
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r814; O=$D/out/runs
mkdir -p $O/K $O/S025 $O/S100 $O/PP025; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_SECULAR_NODES=1 NLEV=2 CREM_NO_DIPOLE_FORCE=1 CREM_PAIR_TORQUE_APSE=1 CREM_PAIR_TORQUE_PS=1 PAIR=electron,positron CREM_PAIR_TORQUE_CONSERVATIVE=1"
N5="u5a,0.5120091953982716 u5b,0.5652161499007791 u5c,0.66666666666666663 u5d,0.8125025652578265 u5e,0.9551918755393795"
echo "start $(date +%T)"
{ for nl in $N5; do IFS=, read nn f <<< "$nl"; echo "S025/F1_ops_${nn}_s8192_5 5 $B CREM_SEPARATION_FLOOR_SCALE=0.025 CREM_MECHANICAL_STEPS_PER_ORBIT=8192 LFRAC=$f PHEN=2"; done
  for h in 4096 2048; do
    for nl in $N5; do IFS=, read nn f <<< "$nl"; echo "S025/F1_ops_${nn}_s${h}_5 5 $B CREM_SEPARATION_FLOOR_SCALE=0.025 CREM_MECHANICAL_STEPS_PER_ORBIT=$h LFRAC=$f PHEN=2"; done
    for nl in $N5; do IFS=, read nn f <<< "$nl"; echo "S100/F1_ops_${nn}_s${h}_5 5 $B CREM_SEPARATION_FLOOR_SCALE=0.1 CREM_MECHANICAL_STEPS_PER_ORBIT=$h LFRAC=$f PHEN=2"; done
  done
  for nl in $N5; do IFS=, read nn f <<< "$nl"; echo "PP025/F0_pps_${nn}_s2048_5 5 $B CREM_SEPARATION_FLOOR_SCALE=0.025 CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=$f PHEN=1"; done
  echo "K/F1_ops_u5a_s2048_5 5 $B CREM_SEPARATION_FLOOR_SCALE=0.05 CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=0.5120091953982716 PHEN=2"
} | xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p814 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

#!/bin/bash
# Audit 812: apsidal precession of the Ps orbit integrated over L (u = 1/L Gauss nodes, n = 2), method of 808; probe p812 (= p808).
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r812; O=$D/out/runs
mkdir -p $O/K $O/P3 $O/P4 $O/PP; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_SECULAR_NODES=1 NLEV=2 CREM_NO_DIPOLE_FORCE=1 CREM_PAIR_TORQUE_APSE=1"
N3="u3a,0.52985793589488495 u3b,0.66666666666666663 u3c,0.89871349267654366"
N4="u4a,0.51798223077348682 u4b,0.59880579377581356 u4c,0.75187434103678785 u4d,0.93507595217092121"
PS="CREM_PAIR_TORQUE_PS=1 PAIR=electron,positron"
echo "start $(date +%T)"
{ for h in 4096 2048 1024; do for nl in $N3; do IFS=, read nn f <<< "$nl"; echo "P3/F1_ops_${nn}_s${h}_5 5 $B $PS CREM_MECHANICAL_STEPS_PER_ORBIT=$h LFRAC=$f PHEN=2"; done
    for nl in $N4; do IFS=, read nn f <<< "$nl"; echo "P4/F1_ops_${nn}_s${h}_5 5 $B $PS CREM_MECHANICAL_STEPS_PER_ORBIT=$h LFRAC=$f PHEN=2"; done
    if [ $h = 2048 ]; then for nl in $N3; do IFS=, read nn f <<< "$nl"; echo "PP/F0_pps_${nn}_s2048_5 5 $B $PS CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=$f PHEN=1"; done
      echo "K/F1_proton_u3b_s2048_5 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=0.66666666666666663 PHEN=2 PAIR=proton,electron"; fi; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p812 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

#!/bin/bash
# Audit 804: control runs of the apsis-detector hysteresis (probe p804): R1/R2 = 803 n4L4 seed 5 (proton has the end-of-window false apocentre), R3 = 799 A n3L1, R4 = 800 B n2L2.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r804; O=$D/out/runs
mkdir -p $O; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1"
B803="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_MECHANICAL_STEPS_PER_ORBIT=512 CREM_SECULAR_NODES=1"
echo "start $(date +%T)"
{ echo "R1_F1_proton_n4L4_5 5 $B803 CREM_PAIR_TORQUE=200 NLEV=4 LFRAC=1 PHEN=2 PAIR=proton,electron"
  echo "R2_F1_antiproton_n4L4_5 5 $B803 CREM_PAIR_TORQUE=200 NLEV=4 LFRAC=1 PHEN=2 PAIR=antiproton,positron"
  echo "R3_proton_n3L1_5 5 $B CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_MECHANICAL_STEPS_PER_ORBIT=512 NLEV=3 LFRAC=0.3333333333333333 PHEN=2 PAIR=proton,electron"
  echo "R4_proton_n2L2_7 7 $B CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_MECHANICAL_STEPS_PER_ORBIT=512 NLEV=2 LFRAC=1 PHEN=2 PAIR=proton,electron"; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p804 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

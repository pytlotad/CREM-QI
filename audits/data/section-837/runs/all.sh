#!/bin/bash
# Audit 837: nature of the 835 residual -- mode ZT of p814 with the start shift split by part (x / v) and particle (1 / 2),
# lambda* per system, and the step 4096 at u5a.  p814 from the tree.
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s837; D=$R/audits/data/section-837; O=$D/out/runs
rm -rf $S; mkdir -p $S $O/K $O/X $O/V $O/E1 $O/E2 $O/LS $O/N
( cd $R && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-814/p814.cpp -o $S/p814 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_SECULAR_NODES=1 NLEV=2 CREM_PAIR_TORQUE_APSE=1 CREM_PAIR_TORQUE_CONSERVATIVE=1"
PS="CREM_PAIR_TORQUE_PS=1 PAIR=electron,positron"; HY="PAIR=proton,electron"; X=CREM_THOMAS_TULCZYJEW_INSTANT=1; T=CREM_THOMAS_TULCZYJEW_START
N5="u5a,0.5120091953982716 u5b,0.5652161499007791 u5c,0.66666666666666663 u5d,0.8125025652578265 u5e,0.9551918755393795"; FC=0.66666666666666663; FA=0.5120091953982716
{ echo "K/D_F1_ops_u5c_s2048 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=2048 $PS $X $T=1 ${T}_MASK=xv12 LFRAC=$FC PHEN=2"
  for nl in $N5; do IFS=, read nn f <<< "$nl"
    for v in "X x12" "V v12" "E1 xv1" "E2 xv2"; do set -- $v; V=$1; Y="$X $T=1 ${T}_MASK=$2"
      echo "$V/D_F1_ops_${nn}_s2048 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=2048 $PS $Y LFRAC=$f PHEN=2"
      echo "$V/D_F0_H_${nn}_s2048 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=2048 $HY $Y LFRAC=$f PHEN=1"
      echo "$V/D_F1_H_${nn}_s2048 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=2048 $HY $Y LFRAC=$f PHEN=2"
    done
    echo "LS/D_F1_ops_${nn}_s2048 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=2048 $PS $X $T=1.0138 LFRAC=$f PHEN=2"
    echo "LS/D_F0_H_${nn}_s2048 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=2048 $HY $X $T=1.0350 LFRAC=$f PHEN=1"
    echo "LS/D_F1_H_${nn}_s2048 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=2048 $HY $X $T=1.1194 LFRAC=$f PHEN=2"
  done
  for sys in "F1_ops PHEN=2 $PS" "F0_H PHEN=1 $HY" "F1_H PHEN=2 $HY"; do set -- $sys; nm=$1; shift
    echo "N/Z_${nm}_u5a_s4096 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=4096 $* LFRAC=$FA"
    echo "N/ZS_${nm}_u5a_s4096 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=4096 $* $X $T=1 LFRAC=$FA"
  done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p814 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

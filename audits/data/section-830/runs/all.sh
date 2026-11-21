#!/bin/bash
# Audit 830: spin arm (charge-moment forces on), retarded (O) vs conservative (Z), 5 Gauss nodes; probe p814 from the tree.  Waits for 829 ALLDONE.
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s830; D=$R/audits/data/section-830; O=$D/out/runs
until grep -q ALLDONE $R/audits/data/section-829/out/runs829.log; do sleep 60; done
echo "829 done $(date +%T)"
rm -rf $S; mkdir -p $S $O/Z $O/O
( cd $R && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-814/p814.cpp -o $S/p814 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_SECULAR_NODES=1 NLEV=2 CREM_PAIR_TORQUE_APSE=1"
PS="CREM_PAIR_TORQUE_PS=1 PAIR=electron,positron"; HY="PAIR=proton,electron"
N4="u5a,0.5120091953982716 u5b,0.5652161499007791 u5d,0.8125025652578265 u5e,0.9551918755393795"; FC=0.66666666666666663
{ echo "Z/D_F1_ops_u5c_s2048 5 $B $PS CREM_PAIR_TORQUE_CONSERVATIVE=1 CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=$FC PHEN=2"
  for m in Z O; do CV=""; [ $m = Z ] && CV=CREM_PAIR_TORQUE_CONSERVATIVE=1
    for nl in $N4; do IFS=, read nn f <<< "$nl"
      echo "$m/D_F1_ops_${nn}_s2048 5 $B $PS $CV CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=$f PHEN=2"
      echo "$m/D_F0_pps_${nn}_s2048 5 $B $PS $CV CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=$f PHEN=1"
      echo "$m/D_F1_H_${nn}_s2048 5 $B $HY $CV CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=$f PHEN=2"
      echo "$m/D_F0_H_${nn}_s2048 5 $B $HY $CV CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=$f PHEN=1"
    done
    echo "$m/D_F1_ops_u5c_s4096 5 $B $PS $CV CREM_MECHANICAL_STEPS_PER_ORBIT=4096 LFRAC=$FC PHEN=2"
    echo "$m/D_F0_H_u5c_s4096 5 $B $HY $CV CREM_MECHANICAL_STEPS_PER_ORBIT=4096 LFRAC=$FC PHEN=1"
  done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p814 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

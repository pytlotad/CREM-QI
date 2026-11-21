#!/bin/bash
# Audit 828: retarded (O) vs conservative Coulomb + Darwin (Z) apsidal rate n = 2; probe p814 from the tree.  Waits for 827 to finish.
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s828; D=$R/audits/data/section-828; O=$D/out/runs
while pgrep -f 'section-827/runs/all.sh|bash runs/all.sh' >/dev/null && ! grep -q ALLDONE $R/audits/data/section-827/out/runs827.log; do sleep 60; done
echo "827 done $(date +%T)"
rm -rf $S; mkdir -p $S $O/Z $O/O
( cd $R && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-814/p814.cpp -o $S/p814 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_SECULAR_NODES=1 NLEV=2 CREM_PAIR_TORQUE_APSE=1"
ND=CREM_NO_DIPOLE_FORCE=1; PS="CREM_PAIR_TORQUE_PS=1 PAIR=electron,positron"; HY="PAIR=proton,electron"
N5="u5a,0.5120091953982716 u5b,0.5652161499007791 u5c,0.66666666666666663 u5d,0.8125025652578265 u5e,0.9551918755393795"
{ for m in Z O; do CV=""; [ $m = Z ] && CV=CREM_PAIR_TORQUE_CONSERVATIVE=1
    for nl in $N5; do IFS=, read nn f <<< "$nl"
      for h in 2048 4096; do echo "$m/F1_ops_${nn}_s${h} 5 $B $ND $PS $CV CREM_MECHANICAL_STEPS_PER_ORBIT=$h LFRAC=$f PHEN=2"; done
      echo "$m/F0_pps_${nn}_s2048 5 $B $ND $PS $CV CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=$f PHEN=1"
      echo "$m/F1_H_${nn}_s2048 5 $B $ND $HY $CV CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=$f PHEN=2"
      echo "$m/F0_H_${nn}_s2048 5 $B $ND $HY $CV CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=$f PHEN=1"
    done
    f=0.66666666666666663
    echo "$m/D_F1_ops_u5c_s2048 5 $B $PS $CV CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=$f PHEN=2"
    echo "$m/D_F0_pps_u5c_s2048 5 $B $PS $CV CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=$f PHEN=1"
    echo "$m/D_F1_H_u5c_s2048 5 $B $HY $CV CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=$f PHEN=2"
    echo "$m/D_F0_H_u5c_s2048 5 $B $HY $CV CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=$f PHEN=1"
  done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p814 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

#!/bin/bash
# Audit 867: proper dipole + BMT factor (841) in mode O through CREM_RETARDED_SPIN_DIRAC=1 (OD) and in the conservative Tulczyjew sum
# (ZD = CONSERVATIVE + THOMAS_TULCZYJEW_INSTANT + START + SPIN_ORBIT_DIRAC) x {o-Ps, p-Ps, H F0, H F1} x 5 nodes; controls K_O, K_Z.
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s867; D=$R/audits/data/section-867; O=$D/out/runs
rm -rf $S; mkdir -p $S $O/K $O/OD $O/ZD
( cd $R && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-814/p814.cpp -o $S/p814 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_SECULAR_NODES=1 NLEV=2 CREM_PAIR_TORQUE_APSE=1 CREM_MECHANICAL_STEPS_PER_ORBIT=2048"
PS="CREM_PAIR_TORQUE_PS=1 PAIR=electron,positron"; HY="PAIR=proton,electron"
XO="CREM_RETARDED_SPIN_DIRAC=1"; XZ="CREM_PAIR_TORQUE_CONSERVATIVE=1 CREM_THOMAS_TULCZYJEW_INSTANT=1 CREM_THOMAS_TULCZYJEW_START=1 CREM_SPIN_ORBIT_DIRAC=1"
E=0.9551918755393795
N5="u5a,0.5120091953982716 u5b,0.5652161499007791 u5c,0.66666666666666663 u5d,0.8125025652578265 u5e,0.9551918755393795"
{ echo "K/K_O 5 $B $PS LFRAC=$E PHEN=2"
  echo "K/K_Z 5 $B $HY CREM_PAIR_TORQUE_CONSERVATIVE=1 LFRAC=$E PHEN=1"
  for nl in $N5; do IFS=, read nn f <<< "$nl"
    for M in OD ZD; do [ $M = ZD ] && X=$XZ || X=$XO
      echo "$M/D_F1_ops_${nn}_s2048 5 $B $PS $X LFRAC=$f PHEN=2"
      echo "$M/D_F0_pps_${nn}_s2048 5 $B $PS $X LFRAC=$f PHEN=1"
      echo "$M/D_F0_H_${nn}_s2048 5 $B $HY $X LFRAC=$f PHEN=1"
      echo "$M/D_F1_H_${nn}_s2048 5 $B $HY $X LFRAC=$f PHEN=2"
    done
  done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p814 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

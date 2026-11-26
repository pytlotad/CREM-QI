#!/bin/bash
# Audit 865: noise floor of mode O in Ps -- start scaled by lambda (CREM_THOMAS_TULCZYJEW_START), no code change:
# O_p5 / O_m5 / O_p6 (mode O, lambda 1.00001 / 0.99999 / 1.000001), Z_p5 (conservative, 1.00001) x {o-Ps, p-Ps} x 5 nodes; control K vs 836 OS.
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s865; D=$R/audits/data/section-865; O=$D/out/runs
rm -rf $S; mkdir -p $S $O/K $O/O_p5 $O/O_m5 $O/O_p6 $O/Z_p5
( cd $R && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-814/p814.cpp -o $S/p814 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_SECULAR_NODES=1 NLEV=2 CREM_PAIR_TORQUE_APSE=1 CREM_MECHANICAL_STEPS_PER_ORBIT=2048"
PS="CREM_PAIR_TORQUE_PS=1 PAIR=electron,positron"; HY="PAIR=proton,electron"; C=CREM_PAIR_TORQUE_CONSERVATIVE=1; T2=CREM_THOMAS_TULCZYJEW_START_ORDER=2
E=0.9551918755393795
N5="u5a,0.5120091953982716 u5b,0.5652161499007791 u5c,0.66666666666666663 u5d,0.8125025652578265 u5e,0.9551918755393795"
N5="u5a,0.5120091953982716 u5b,0.5652161499007791 u5c,0.66666666666666663 u5d,0.8125025652578265 u5e,0.9551918755393795"
{ echo "K/K_O 5 $B $PS CREM_THOMAS_TULCZYJEW_START=1 LFRAC=0.5120091953982716 PHEN=2"
  for nl in $N5; do IFS=, read nn f <<< "$nl"
    for V in O_p5,,1.00001 O_m5,,0.99999 O_p6,,1.000001 Z_p5,$C,1.00001; do IFS=, read M X L <<< "$V"
      echo "$M/D_F1_ops_${nn}_s2048 5 $B $PS $X CREM_THOMAS_TULCZYJEW_START=$L LFRAC=$f PHEN=2"
      echo "$M/D_F0_pps_${nn}_s2048 5 $B $PS $X CREM_THOMAS_TULCZYJEW_START=$L LFRAC=$f PHEN=1"
    done
  done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p814 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

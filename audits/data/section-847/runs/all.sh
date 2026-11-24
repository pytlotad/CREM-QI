#!/bin/bash
# Audit 847: controls of the ZDC default (new code): N (defaults = 846), L (CREM_SPIN_FORM_LEGACY = 830), T (Tulczyjew research = 835 / 836).
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s847; D=$R/audits/data/section-847; O=$D/out/runs
rm -rf $S; mkdir -p $S $O
( cd $R && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-814/p814.cpp -o $S/p814 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_SECULAR_NODES=1 NLEV=2 CREM_PAIR_TORQUE_APSE=1 CREM_MECHANICAL_STEPS_PER_ORBIT=2048"
PS="CREM_PAIR_TORQUE_PS=1 PAIR=electron,positron"; HY="PAIR=proton,electron"; C=CREM_PAIR_TORQUE_CONSERVATIVE=1; LG=CREM_SPIN_FORM_LEGACY=1
E=0.9551918755393795; FC=0.66666666666666663
{ echo "N_O 5 $B $PS LFRAC=$E PHEN=2"
  echo "L_O 5 $B $PS $LG LFRAC=$E PHEN=2"
  echo "T_OS 5 $B $PS CREM_THOMAS_TULCZYJEW_START=1 LFRAC=$E PHEN=2"
  echo "N_Z 5 $B $PS $C LFRAC=$E PHEN=2"
  echo "N_ZH 5 $B $HY $C LFRAC=$E PHEN=1"
  echo "L_Z 5 $B $PS $C $LG LFRAC=$FC PHEN=2"
  echo "T_ZS 5 $B $PS $C CREM_THOMAS_TULCZYJEW_INSTANT=1 CREM_THOMAS_TULCZYJEW_START=1 LFRAC=$E PHEN=2"; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p814 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

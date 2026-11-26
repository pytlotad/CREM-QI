#!/bin/bash
# Audit 863: controls of the O-form default of p814 (new code): N_O = 836 OS, N_Z = 835 ZS (o-Ps, p-Ps, H F0, H F1), V (ZDC_LEGACY) = 846,
# L (LEGACY) = 830, E_O (explicit 846) = 846 OZ, E_ZT (explicit as 862 ZDCT) = 862; production sets 2 and 5Ps = 861 archive.
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s863; D=$R/audits/data/section-863; O=$D/out/runs
rm -rf $S; mkdir -p $S $O
( cd $R && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-814/p814.cpp -o $S/p814 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
cd $R && rsync -a --exclude .git --exclude audits ./ $S/prod/ && ( cd $S/prod && make -j3 crem_interactions >/dev/null 2>&1 ) || { echo "BLAD produkcja"; exit 1; }
echo "build $(date +%T)"; export S O
for t in "2 --mode statistical --phenomenon 2 --runs 16 --seed 466" "5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551"; do set -- $t; tag=$1; shift
  mkdir -p $S/reg/$tag; ( ulimit -v 2500000; cd $S/reg/$tag; $S/prod/crem_interactions "$@" > out.txt 2> err.txt; echo $? > rc.txt ); echo "prod $tag $(date +%T)"; done
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_SECULAR_NODES=1 NLEV=2 CREM_PAIR_TORQUE_APSE=1 CREM_MECHANICAL_STEPS_PER_ORBIT=2048"
PS="CREM_PAIR_TORQUE_PS=1 PAIR=electron,positron"; HY="PAIR=proton,electron"; C=CREM_PAIR_TORQUE_CONSERVATIVE=1
E=0.9551918755393795; FC=0.66666666666666663
{ for M in O Z; do [ $M = Z ] && X=$C || X=""
    echo "N_${M}_ops 5 $B $PS $X LFRAC=$E PHEN=2"
    echo "N_${M}_pps 5 $B $PS $X LFRAC=$E PHEN=1"
    echo "N_${M}_H0 5 $B $HY $X LFRAC=$E PHEN=1"
    echo "N_${M}_H1 5 $B $HY $X LFRAC=$E PHEN=2"; done
  echo "V_O 5 $B $PS CREM_SPIN_FORM_ZDC_LEGACY=1 LFRAC=$E PHEN=2"
  echo "V_Z 5 $B $PS $C CREM_SPIN_FORM_ZDC_LEGACY=1 LFRAC=$E PHEN=2"
  echo "L_O 5 $B $PS CREM_SPIN_FORM_LEGACY=1 LFRAC=$E PHEN=2"
  echo "L_Z 5 $B $PS $C CREM_SPIN_FORM_LEGACY=1 LFRAC=$FC PHEN=2"
  echo "E_O 5 $B $PS CREM_RETARDED_SPIN_ZDC=1 LFRAC=$E PHEN=2"
  echo "E_ZT 5 $B $PS $C CREM_SPIN_ORBIT_DIRAC=1 CREM_SPIN_ORBIT_CANONICAL=1 CREM_THOMAS_TULCZYJEW_START=1 LFRAC=$E PHEN=2"; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p814 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

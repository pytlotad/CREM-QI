#!/bin/bash
# Audit 843: CREM_SPIN_ORBIT_CANONICAL on top of 841 -- static S (p843) and P (p843b),
# then p814 ZDCA (= 841 ZDA + CANONICAL) x {o-Ps, H F0, H F1} x 5 nodes + K (ZDA H F0 u5c, 841 only).
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s843; D=$R/audits/data/section-843; O=$D/out/runs
rm -rf $S; mkdir -p $S $O/S $O/P $O/K $O/ZDCA
cd $R
for p in audits/data/section-814/p814 audits/data/section-843/p843 audits/data/section-843/p843b; do
  g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules $p.cpp -o $S/$(basename $p) $(root-config --libs) 2>/dev/null || { echo "BLAD kompilacja $p"; exit 1; }
done
echo "build $(date +%T)"; export S O
declare -A LF=([u5a]=0.5120091953982716 [u5b]=0.5652161499007791 [u5c]=0.66666666666666663 [u5d]=0.8125025652578265 [u5e]=0.9551918755393795)
while read sk nd pair rp e a b c x y z; do
  ( ulimit -v 2500000; $S/p843b $pair $rp $e $a $b $c $x $y $z 4096 > $O/P/${sk}_${nd}.txt; $S/p843b $pair $rp $e $a $b $c $x $y $z 8192 > $O/P/${sk}_${nd}_N2.txt )
  if [ $sk != ops ]; then ( ulimit -v 2500000; $S/p843 ${LF[$nd]} $x $y $z 4096 65536 R > $O/S/${sk}_${nd}.txt; $S/p843 ${LF[$nd]} $x $y $z 4096 131072 R > $O/S/${sk}_${nd}_h2.txt ); fi
done < $D/out/params839.txt
echo "static $(date +%T)"
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_SECULAR_NODES=1 NLEV=2 CREM_PAIR_TORQUE_APSE=1 CREM_PAIR_TORQUE_CONSERVATIVE=1 CREM_MECHANICAL_STEPS_PER_ORBIT=2048 CREM_PAIR_TORQUE_START_APO=1"
PS="CREM_PAIR_TORQUE_PS=1 PAIR=electron,positron"; HY="PAIR=proton,electron"; X="CREM_SPIN_ORBIT_DIRAC=1 CREM_SPIN_ORBIT_CANONICAL=1"
{ echo "K/D_F0_H_u5c_s2048 5 $B $HY CREM_SPIN_ORBIT_DIRAC=1 LFRAC=0.66666666666666663 PHEN=1"
  for nd in u5a u5b u5c u5d u5e; do f=${LF[$nd]}
    echo "ZDCA/D_F1_ops_${nd}_s2048 5 $B $PS $X LFRAC=$f PHEN=2"
    echo "ZDCA/D_F0_H_${nd}_s2048 5 $B $HY $X LFRAC=$f PHEN=1"
    echo "ZDCA/D_F1_H_${nd}_s2048 5 $B $HY $X LFRAC=$f PHEN=2"
  done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p814 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

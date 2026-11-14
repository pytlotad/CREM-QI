#!/bin/bash
# Audit 789: 44 clock runs (part A) + 6 torque runs (part B); probe = ../section-710/p710.cpp built from HEAD 2801873e as r789/p789.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r789; O=$D/out/runs
mkdir -p $O/P $O/K/ret $O/K/ins $O/X3/ret $O/X3/ins $O/XH/ret $O/XH/ins $O/B; export S O
CK="MODE=est QSPIN=1 BUDGET=5 CREM_PAIR_CLOCK=20 CREM_MECHANICAL_STEPS_PER_ORBIT=8192"
INS="CREM_PAIR_CLOCK_INSTANT=1 CREM_MECHANICAL_HAZARD_SCALE=0"
HY="CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 PHEN=2"
echo "start $(date +%T)"
{ for s in 2 3; do for ph in 1 2; do for e in 0 1e-3 -1e-3 2e-3 -2e-3; do
    v="PHEN=$ph"; [ "$e" != 0 ] && v="PHEN=$ph CREM_PAIR_CLOCK_VSCALE=$e"; echo "P/ps_A${ph}_${s}_e$e $s $CK NLEV=2 LFRAC=0.75 $v"; done; done; done
  for sm in ret ins; do x="QSPIN=1"; [ $sm = ins ] && x="$INS"
    for s in 0 1; do for ph in 1 2; do
      echo "K/$sm/ps_A${ph}_${s}_e0 $s $CK NLEV=2 LFRAC=0.75 PHEN=$ph $x"
      echo "X3/$sm/ps_A${ph}_${s}_e0 $s $CK NLEV=3 LFRAC=0.8333333333333334 PHEN=$ph $x"; done; done
    for s in 2 3; do echo "XH/$sm/proton_${s}_e0 $s $CK NLEV=2 LFRAC=0.75 PAIR=proton,electron $HY $x"
      echo "XH/$sm/antiproton_${s}_e0 $s $CK NLEV=2 LFRAC=0.75 PAIR=antiproton,positron $HY $x"; done; done
  for s in 0 2 3; do for t in "proton PAIR=proton,electron" "antiproton PAIR=antiproton,positron"; do set -- $t
    echo "B/${1}_$s $s MODE=est QSPIN=1 BUDGET=60 NLEV=2 LFRAC=0.75 CREM_PAIR_TORQUE=200 $2 $HY"; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p789 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

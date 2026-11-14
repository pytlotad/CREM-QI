#!/bin/bash
# Audit 789, POPRAWKA 2: the 28 eps-0 clock runs again with CREM_PAIR_CLOCK_VSCALE=0 (prints Lq / LTh; dynamics unchanged).
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r789; O=$D/out/runs; export S O
CK="MODE=est QSPIN=1 BUDGET=5 CREM_PAIR_CLOCK=20 CREM_MECHANICAL_STEPS_PER_ORBIT=8192 CREM_PAIR_CLOCK_VSCALE=0"
INS="CREM_PAIR_CLOCK_INSTANT=1 CREM_MECHANICAL_HAZARD_SCALE=0"
HY="CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 PHEN=2"
echo "start $(date +%T)"
{ for s in 2 3; do for ph in 1 2; do echo "P/ps_A${ph}_${s}_e0 $s $CK NLEV=2 LFRAC=0.75 PHEN=$ph"; done; done
  for sm in ret ins; do x="QSPIN=1"; [ $sm = ins ] && x="$INS"
    for s in 0 1; do for ph in 1 2; do
      echo "K/$sm/ps_A${ph}_${s}_e0 $s $CK NLEV=2 LFRAC=0.75 PHEN=$ph $x"
      echo "X3/$sm/ps_A${ph}_${s}_e0 $s $CK NLEV=3 LFRAC=0.8333333333333334 PHEN=$ph $x"; done; done
    for s in 2 3; do echo "XH/$sm/proton_${s}_e0 $s $CK NLEV=2 LFRAC=0.75 PAIR=proton,electron $HY $x"
      echo "XH/$sm/antiproton_${s}_e0 $s $CK NLEV=2 LFRAC=0.75 PAIR=antiproton,positron $HY $x"; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p789 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

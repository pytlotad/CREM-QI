#!/bin/bash
# Audit 792: K (1 control), T (2 chunk traces, k = 1), S (8 runs, CREM_SPIN_RATE_SECTOR=charge), D (2 runs, CREM_NO_DIPOLE_FORCE=1).
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r792; O=$D/out/runs
mkdir -p $O/K $O/T $O/S $O/D; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 PHEN=2 NLEV=3 CREM_PAIR_TORQUE_ANYL=1 CREM_MECHANICAL_STEPS_PER_ORBIT=512"
L1="LFRAC=0.3333333333333333"; L15="LFRAC=0.5"
H="PAIR=proton,electron"; A="PAIR=antiproton,positron"
echo "start $(date +%T)"
{ echo "T/proton_n3L1_5 5 $B $L1 $H CREM_PAIR_TORQUE_CHUNK=1"
  echo "T/antiproton_n3L1_5 5 $B $L1 $A CREM_PAIR_TORQUE_CHUNK=1"
  echo "K/proton_n3L1_5 5 $B $L1 $H"
  for s in 4 5 6; do echo "S/proton_n3L1_$s $s $B $L1 $H CREM_SPIN_RATE_SECTOR=charge"; echo "S/antiproton_n3L1_$s $s $B $L1 $A CREM_SPIN_RATE_SECTOR=charge"; done
  echo "S/proton_n3L1.5_5 5 $B $L15 $H CREM_SPIN_RATE_SECTOR=charge"; echo "S/antiproton_n3L1.5_5 5 $B $L15 $A CREM_SPIN_RATE_SECTOR=charge"
  echo "D/proton_n3L1_5 5 $B $L1 $H CREM_NO_DIPOLE_FORCE=1"; echo "D/antiproton_n3L1_5 5 $B $L1 $A CREM_NO_DIPOLE_FORCE=1"; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p792 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

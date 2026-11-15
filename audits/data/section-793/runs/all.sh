#!/bin/bash
# Audit 793: E (6 energy traces), Z (2 convergence runs), G (48 rule runs with CREM_NO_DIPOLE_FORCE=1).
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r793; O=$D/out/runs
mkdir -p $O/E $O/Z $O/G; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 PHEN=2 CREM_PAIR_TORQUE_ANYL=1"
TR="NLEV=3 LFRAC=0.3333333333333333 CREM_PAIR_TORQUE_CHUNK=1 CREM_PAIR_TORQUE_ENERGY=1"
echo "start $(date +%T)"
{ for t in "proton PAIR=proton,electron" "antiproton PAIR=antiproton,positron"; do set -- $t
    echo "E/${1}_E512 5 $B $TR $2 CREM_MECHANICAL_STEPS_PER_ORBIT=512"
    echo "E/${1}_E0 5 $B $TR $2 CREM_MECHANICAL_STEPS_PER_ORBIT=512 CREM_NO_DIPOLE_FORCE=1"
    echo "E/${1}_E2048 5 $B $TR $2 CREM_MECHANICAL_STEPS_PER_ORBIT=2048"; done
  echo "Z/proton_n3L1_4 4 $B NLEV=3 LFRAC=0.3333333333333333 PAIR=proton,electron CREM_NO_DIPOLE_FORCE=1 CREM_MECHANICAL_STEPS_PER_ORBIT=2048"
  echo "Z/proton_n2L1_4 4 $B NLEV=2 LFRAC=0.5 PAIR=proton,electron CREM_NO_DIPOLE_FORCE=1 CREM_MECHANICAL_STEPS_PER_ORBIT=2048"
  for nl in "3 L1 0.3333333333333333" "3 L1.5 0.5" "3 L2 0.6666666666666666" "3 L2.5 0.8333333333333334" "3 L3 1" "2 L1 0.5" "2 L1.5 0.75" "2 L2 1"; do
    for s in 4 5 6; do for t in "proton PAIR=proton,electron" "antiproton PAIR=antiproton,positron"; do set -- $t $nl
      echo "G/${1}_n$3$4_$s $s $B NLEV=$3 LFRAC=$5 $2 CREM_NO_DIPOLE_FORCE=1 CREM_MECHANICAL_STEPS_PER_ORBIT=512"; done; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p793 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

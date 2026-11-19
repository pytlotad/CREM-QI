#!/bin/bash
# Audit 808: apsidal precession integrated over L (u = 1/L Gauss nodes, n = 2) + anomaly from BMT (g vs 2 at L = 1.5); probe p808.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r808; O=$D/out/runs
mkdir -p $O/O1 $O/O2 $O/O3 $O/O4 $O/A1 $O/A2; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_SECULAR_NODES=1 NLEV=2 CREM_NO_DIPOLE_FORCE=1 CREM_PAIR_TORQUE_APSE=1"
N3="u3a,0.52985793589488495 u3b,0.66666666666666663 u3c,0.89871349267654366"
N4="u4a,0.51798223077348682 u4b,0.59880579377581356 u4c,0.75187434103678785 u4d,0.93507595217092121"
declare -A TYP=([proton]="PAIR=proton,electron" [antiproton]="PAIR=antiproton,positron")
job() { echo "$1/F$2_$3_$4_s$5_$6 $6 $B CREM_MECHANICAL_STEPS_PER_ORBIT=$5 LFRAC=$7 PHEN=$(($2 + 1)) ${TYP[$3]}${8:+ $8}"; }
echo "start $(date +%T)"
{ for F in 1 0; do for t in proton antiproton; do for g in "gE" "g2 CREM_PAIR_TORQUE_BMT_G=2"; do set -- $g
      job A1 $F $t $1 2048 5 0.75 "${2:-}"; done; done; done
  for F in 1 0; do for t in proton antiproton; do for nl in $N3; do IFS=, read nn f <<< "$nl"; job O1 $F $t $nn 2048 5 $f; done; done; done
  for nl in $N4; do IFS=, read nn f <<< "$nl"; job O3 1 proton $nn 2048 5 $f; done
  job O4 1 proton u3b 2048 6 0.66666666666666663
  for h in 1024 512; do for nl in $N3; do IFS=, read nn f <<< "$nl"; job O2 1 proton $nn $h 5 $f; done; done
  for g in "gE" "g2 CREM_PAIR_TORQUE_BMT_G=2"; do set -- $g; job A2 1 proton $1 1024 5 0.75 "${2:-}"; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p808 $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

#!/bin/bash
# Audit 868 (prototype quasi-3D, E1): guide868.cpp; arms C (no Q3D, = 773 G_01 first 300 records), E (Q3D, beta 0, norm/levels),
# Q (Q3D, beta 0.2, CAUST = CAUSTT = -pi/2, 773/01-40, 12000 periods), QT (CAUSTT 0, 773/41-60, 6000), QN (beta -0.05, 773/01-40, 12000).
# Waits for audit 867 to finish (audits not in parallel).  3 processes, ulimit -v 2500000 each, LC_ALL=C.
set -u; export LC_ALL=C; ulimit -v 4000000
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s868; B=$S/guide868
until grep -q "ALLDONE\|BLAD" $D/../section-867/out/runs867.log; do sleep 120; done
echo "867 done, start $(date +%T)"
mkdir -p $S; O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/guide868.cpp || { echo BLAD kompilacja; exit 1; }
export B O
ST=$D/../section-773/starts773.txt; PI2=-1.5707963267948966
{ echo "C 0 0.2 0.001 $PI2 0 300 01 $(awk '$1=="01"{print $2}' $ST)"
  echo "E 1 0 0.001 $PI2 $PI2 1 01 $(awk '$1=="01"{print $2}' $ST)"
  awk -v p=$PI2 '$1<=40{print "Q 1 0.2 0.001 "p" "p" 12000 "$1" "$2}' $ST
  awk -v p=$PI2 '$1>40{print "QT 1 0.2 0.001 "p" 0 6000 "$1" "$2}' $ST
  awk -v p=$PI2 '$1<=40{print "QN 1 -0.05 0.001 "p" "p" 12000 "$1" "$2}' $ST; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; Q3D=$1 BETA=$2 GT=$3 CAUST=$4 CAUSTT=$5 A=1 INIT="$8" $B $6 0 500 A > $O/${0}_$7.txt && echo "done ${0}_$7 $(date +%T)" || echo "FAIL ${0}_$7 $(date +%T)"'
echo "ALLDONE $(date +%T)"

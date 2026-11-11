#!/bin/bash
# Audit 780 runs (scale test): K control, A step convergence (pure modes), B quality GT 0.0001, C ratio beta 0.0125 / GT 0.000075.  Model ../section-772/guide772.cpp unchanged.
# 3 parallel (xargs -P 3), ulimit -v 4 GB each, LC_ALL=C.  Longest jobs first.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s780/guide772
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/../section-772/guide772.cpp || { echo BLAD kompilacja; exit 1; }
export B O; ST=$D/../section-773/starts773.txt
( ulimit -v 4000000; BETA=0.2 GT=0.0003 CAUST=-1.5707963267948966 A=1 INIT="$(awk '$1=="01"{print $2}' $ST)" $B 300 0 500 A > $O/K_ctrl.txt ) && echo "done K $(date +%T)"
{ awk '$1<=20{print "C_"$1" 0.0125 0.000075 2000 40000 "$2}' $ST
  awk '$1<=20{print "B_"$1" 0.2 0.0001 2000 30000 "$2}' $ST
  for s in 8000 4000 2000; do for t in "0.05 1 0" "0.05 2 0" "0.2 2 0" "0.05 0 0" "0.05 0 2" "0.2 1 0" "0.2 0 2"; do set -- $t
    echo "A_b$1_$2_$3_s$s $1 0.0003 $s 4000 $2:$3:1:0"; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 4000000; BETA=$1 GT=$2 STEPS=$3 CAUST=-1.5707963267948966 A=1 INIT="$5" $B $4 0 500 A > $O/$0.txt && echo "done $0 $(date +%T)"'
echo "ALLDONE $(date +%T)"

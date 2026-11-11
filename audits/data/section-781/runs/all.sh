#!/bin/bash
# Audit 781 runs: K controls (done first, also here), P pure modes RAD 1/10 (10000 periods), R1/R10 starts 773/01-20 (20000 periods).  beta 0.2, GT 0.0003, CAUST -pi/2.
# 3 parallel (xargs -P 3), ulimit -v 4 GB each, LC_ALL=C.  Longest first.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s781/guide781
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/guide781.cpp || { echo BLAD kompilacja; exit 1; }
export B O; ST=$D/../section-773/starts773.txt; I="$(awk '$1=="01"{print $2}' $ST)"
( ulimit -v 4000000; BETA=0.2 GT=0.0003 CAUST=-1.5707963267948966 A=1 INIT="$I" $B 300 0 500 A > $O/K_norad.txt
  BETA=0.2 GT=0.0003 CAUST=-1.5707963267948966 A=1 RAD=0 INIT="$I" $B 300 0 500 A > $O/K_rad0.txt ) && echo "done K $(date +%T)"
{ for r in 1 10; do awk -v r=$r '$1<=20{print "R"r"_"$1" "r" 20000 "$2}' $ST; done
  for r in 1 10; do for md in 0:0 0:1 0:2 1:1 1:-2 2:0; do echo "P${r}_${md/:/_} $r 10000 $md:1:0"; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 4000000; BETA=0.2 GT=0.0003 CAUST=-1.5707963267948966 A=1 RAD=$1 INIT="$3" $B $2 0 500 A > $O/$0.txt && echo "done $0 $(date +%T)"'
echo "ALLDONE $(date +%T)"

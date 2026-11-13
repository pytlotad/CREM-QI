#!/bin/bash
# Audit 784 runs: E controls (no NORM), N pure modes at P_fix 1e4 / 1e6 / 1e8 (30000 periods, WDUMP), R6 / R8 starts 773/01-20 (15000 periods).
# beta 0.2, CAUST -pi/2, A 1, GT 0, GREF 0.0003, DIP 30.  3 parallel (xargs -P 3), ulimit -v 4 GB each, LC_ALL=C.  Longest first.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s784/guide784
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/guide784.cpp || { echo BLAD kompilacja; exit 1; }
export B O; ST=$D/../section-773/starts773.txt; I="$(awk '$1=="01"{print $2}' $ST)"
echo "start $(date +%T)"
( ulimit -v 4000000; export BETA=0.2 CAUST=-1.5707963267948966
  GT=0.0003 A=1 INIT="$I" $B 300 0 500 A > $O/E_nonorm.txt
  GT=0 GREF=0.0003 DIP=30 A=1 WDUMP=1 INIT="0:1:1:0" $B 300 0 500 A > $O/E_D.txt ) && echo "done E $(date +%T)"
{ for p in 1e4 1e6 1e8; do for md in 0:1 0:2 1:1 1:-2 2:0 0:0; do echo "N${p}_${md/:/_} $p 30000 1 $md:1:0"; done; done
  for p in 1e6 1e8; do awk -v p=$p '$1<=20{print "R"substr(p,3)"_"$1" "p" 15000 0 "$2}' $ST; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 4000000; export BETA=0.2 CAUST=-1.5707963267948966 GT=0 GREF=0.0003 DIP=30 A=1 NORM=$1 WDUMP=$3 INIT="$4"
  $B $2 0 $1 A > $O/$0.txt && echo "done $0 $(date +%T)"'
echo "ALLDONE $(date +%T)"

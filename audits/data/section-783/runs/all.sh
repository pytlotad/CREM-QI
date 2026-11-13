#!/bin/bash
# Audit 783 runs: E controls, K code check (A = 0, GW 10), P pure modes in arms WD / W / D at s = 30 (30000 periods, WDUMP), R (WD, s = 1) starts 773/01-20 (20000 periods).
# beta 0.2, CAUST -pi/2, GREF 0.0003.  3 parallel (xargs -P 3), ulimit -v 4 GB each, LC_ALL=C.  Longest first.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s783/guide783
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/guide783.cpp || { echo BLAD kompilacja; exit 1; }
export B O; ST=$D/../section-773/starts773.txt; I="$(awk '$1=="01"{print $2}' $ST)"
echo "start $(date +%T)"
( ulimit -v 4000000; export BETA=0.2 CAUST=-1.5707963267948966
  GT=0.0003 A=1 INIT="$I" $B 300 0 500 A > $O/E_nogw.txt
  GT=0.0003 A=1 GW=0 INIT="$I" $B 300 0 500 A > $O/E_gw0.txt
  GT=0.0003 A=0 GW=10 WDUMP=1 INIT="0:1:0.99:0;0:0:0.01:0" $B 3000 0 500 A > $O/K_01_00.txt ) && echo "done E K $(date +%T)"
{ for arm in WD W D; do for md in 0:1 0:2 1:1 1:-2 2:0 0:0; do echo "P${arm}_${md/:/_} $arm 30000 1 $md:1:0 30"; done; done
  awk '$1<=20{print "R_"$1" WD 20000 0 "$2" 1"}' $ST; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 4000000; export BETA=0.2 CAUST=-1.5707963267948966 GREF=0.0003 A=1 WDUMP=$3 INIT="$4"
  case $1 in WD) export GT=0.0003 GW=$5 DIP=$5;; W) export GT=0.0003 GW=$5;; D) export GT=0 DIP=$5;; esac
  $B $2 0 500 A > $O/$0.txt && echo "done $0 $(date +%T)"'
echo "ALLDONE $(date +%T)"

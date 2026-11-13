#!/bin/bash
# Audit 785 runs (engine unchanged: ../section-784/guide784.cpp without NORM): E controls, P pure modes at DIP 30/100/300 (15000 periods, WDUMP), R starts 773/01-20 at DIP 30/100/300 (10000 periods).
# beta 0.2, GT 0.0003, CAUST -pi/2, A 1, P0 500.  3 parallel (xargs -P 3), ulimit -v 4 GB each, LC_ALL=C.  Longest first.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s785/guide
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/../section-784/guide784.cpp || { echo BLAD kompilacja; exit 1; }
export B O; ST=$D/../section-773/starts773.txt; I="$(awk '$1=="01"{print $2}' $ST)"
echo "start $(date +%T)"
( ulimit -v 4000000; export BETA=0.2 CAUST=-1.5707963267948966 GT=0.0003 A=1
  INIT="$I" $B 300 0 500 A > $O/E_nodip.txt
  DIP=10 WDUMP=1 INIT="0:1:1:0" $B 300 0 500 A > $O/E10.txt ) && echo "done E $(date +%T)"
{ for d in 30 100 300; do for md in 0:1 0:2 1:1 1:-2 2:0 0:0; do echo "P${d}_${md/:/_} $d 15000 1 $md:1:0"; done; done
  for d in 30 100 300; do awk -v d=$d '$1<=20{print "R"d"_"$1" "d" 10000 0 "$2}' $ST; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 4000000; export BETA=0.2 CAUST=-1.5707963267948966 GT=0.0003 A=1 DIP=$1 WDUMP=$3 INIT="$4"
  $B $2 0 500 A > $O/$0.txt && echo "done $0 $(date +%T)"'
echo "ALLDONE $(date +%T)"

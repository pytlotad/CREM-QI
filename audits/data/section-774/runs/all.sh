#!/bin/bash
# Audit 774 runs: 3 parallel processes (xargs -P 3), each ulimit -v 4 GB.  LC_ALL=C.  Model: ../section-772/guide772.cpp unchanged; starts ../section-773/starts773.txt.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s774/guide772}
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/../section-772/guide772.cpp || { echo BLAD kompilacja; exit 1; }
export B O
{ awk '$1<=40{print "G 0.0003 -1.5707963267948966 20000 "$1" "$2}' $D/../section-773/starts773.txt
  awk '$1>40{print "K 0.0003 0 10000 "$1" "$2}' $D/../section-773/starts773.txt; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 4000000; GT=$1 CAUST=$2 A=1 INIT="$5" $B $3 0 500 A > $O/${0}_$4.txt && echo "done ${0}_$4"'
echo ALLDONE

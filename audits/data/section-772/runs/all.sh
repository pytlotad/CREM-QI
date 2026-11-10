#!/bin/bash
# Audit 772 runs (sequential).  LC_ALL=C.  Memory cap (766).
set -u; export LC_ALL=C; ulimit -v 4000000
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s772/guide772}
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/guide772.cpp || { echo BLAD kompilacja; exit 1; }
GT=0.1 CAUST=-1.5707963267948966 A=1 $B 300 0.3 500 A > $O/C_ctrl.txt; echo C
while read id ini; do GT=0.001 CAUST=-1.5707963267948966 A=1 INIT="$ini" $B 3000 0 500 A > $O/G_$id.txt; done < $D/starts772.txt; echo G
head -20 $D/starts772.txt | while read id ini; do GT=0.001 CAUST=0 A=1 INIT="$ini" $B 3000 0 500 A > $O/K_$id.txt; done; echo K
head -12 $D/starts772.txt | while read id ini; do GT=0.0003 CAUST=-1.5707963267948966 A=1 INIT="$ini" $B 10000 0 500 A > $O/W_$id.txt; done; echo W
echo ALLDONE

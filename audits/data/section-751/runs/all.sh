#!/bin/bash
# Audit 751 runs (sequential).  LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s751/rad751}
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/rad751.cpp || exit 1
for K in 1e5 1e4 1e3; do $B $K 1e-2 15 > $O/K${K}_s1e-2.txt; done
$B 1e4 0 20 > $O/K1e4_s0.txt
S0=$(cat $D/out/s0.txt); $B 1e4 $S0 15 > $O/K1e4_sorbit.txt
mkdir -p $O/ens; for i in $(seq 0 99); do S=$(python3 -c "print(($i+0.5)/100)"); $B 1e5 $S 12 0.02 1000 > $O/ens/s$i.txt; done
echo done

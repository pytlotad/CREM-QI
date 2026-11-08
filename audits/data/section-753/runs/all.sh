#!/bin/bash
# Audit 753 runs (sequential ensembles; OpenMP inside one ensemble).  LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s753/rad753}
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -fopenmp -std=c++17 -o $B $D/rad753.cpp || exit 1
$B 1e5 0 0 0    1 5 > $O/C_2p_nozpf.txt
$B 1e5 0 0 1e-2 1 15 0.05 1500 0.1 200 > $O/D_s1e-2_nozpf.txt
$B 1e5 0 1 0 1000 5 > $O/A_2p_zpf.txt
$B 1e5 1 1 0 1000 5 > $O/B_1s_zpf.txt
echo done

#!/bin/bash
# Audit 749 runs (sequential).  LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s749/lock749}
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -Wall -o $B $D/lock749.cpp || exit 1
JR0=1.5 L0=0.0 GT=0.1 $B 100 10 0.5 1 0 0 > $O/ctrl_noDUMP_Jr1.5_L0.0.txt
DUMP=1 JR0=1.5 L0=0.0 GT=0.1 $B 100 10 0.5 1 0 0 > $O/A_full_Jr1.5_L0.0.txt
MMAX=0 JR0=1.5 L0=0.0 GT=0.1 $B 100 10 0.5 1 0 0 > $O/B_m0_Jr1.5_L0.0.txt
DUMP=1 MMIN=0 JR0=1.5 L0=0.0 GT=0.1 $B 100 10 0.5 1 0 0 > $O/C_mpos_Jr1.5_L0.0.txt
DUMP=1 JR0=0.5 L0=1.0 GT=0.1 $B 100 10 0.5 1 0 0 > $O/D_full_Jr0.5_L1.0.txt
for K in 30 100; do JR0=1.5 L0=0.0 GT=0.1 $B 100 $K 0.5 1 0 0 > $O/E_k${K}_Jr1.5_L0.0.txt; done
echo done

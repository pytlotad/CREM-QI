#!/bin/bash
# Audit 759 runs (sequential): soliton source (SOL=1).  nu T = 10 kappa_eff T, kappa_a T = gamma T / 10 (kappa_eff = nu kappa_a / gamma).  LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s759/lock759}
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/lock759.cpp || exit 1
G=$(python3 -c "import math;N=0.5+math.sqrt(1.4);print(0.1/(2*math.pi*N**3))")
for KE in 3 10 30; do for GT in 0.03 0.1 0.3 1 10; do
  NUT=$(python3 -c "print(10*$KE)"); KAT=$(python3 -c "print($GT/10)")
  for J in 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 1.00; do
    SOL=1 NUT=$NUT KAT=$KAT GAT=$GT JR0=$J L0=1.0 GABS=$G $B 200 0 0.5 1 0 0 > $O/S_ke${KE}_g${GT}_Jr$J.txt
  done; done; done
echo done

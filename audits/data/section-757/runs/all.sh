#!/bin/bash
# Audit 757 runs (sequential).  LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s757/lock757}
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/lock757.cpp || exit 1
G=$(python3 -c "import math;N=0.5+math.sqrt(1.4);print(0.1/(2*math.pi*N**3))")
JR0=0.50 L0=1.0 GABS=$G $B 60 0 0.5 1 0 0 > $O/ctrl744_k0_Jr0.50.txt
for KT in 10 0; do
  for J in $(seq -f %.2f 0.20 0.05 1.30); do JR0=$J L0=1.0 GABS=$G $B 100 $KT 0.5 1 0 0 > $O/J_k${KT}_Jr$J.txt; done
  for L in $(seq -f %.2f 0.60 0.05 1.40); do JR0=0.5 L0=$L GABS=$G $B 100 $KT 0.5 1 0 0 > $O/L_k${KT}_L$L.txt; done
done
echo done

#!/bin/bash
# Audit 760 runs (sequential).  LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=${S:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s760}
mkdir -p $S; O=$D/out/runs; mkdir -p $O
g++ -O3 -march=native -std=c++17 -o $S/lock757 $D/../section-757/lock757.cpp || exit 1
g++ -O3 -march=native -std=c++17 -o $S/lock759 $D/../section-759/lock759.cpp || exit 1
G=$(python3 -c "import math;N=0.5+math.sqrt(1.4);print(0.1/(2*math.pi*N**3))")
for J in $(seq -f %.3f 0.200 0.025 1.400); do
  for KT in 1 3; do JR0=$J L0=1.0 GABS=$G $S/lock757 200 $KT 0.5 1 0 0 > $O/k${KT}_Jr$J.txt; done
  SOL=1 NUT=100 KAT=0.01 GAT=0.1 JR0=$J L0=1.0 GABS=$G $S/lock759 200 0 0.5 1 0 0 > $O/sol_Jr$J.txt
done
echo done

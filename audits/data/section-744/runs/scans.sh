#!/bin/bash
# Audit 744: re-run of the H744b-d scans only (first run: locale seq produced "0,30" -> JR0 parsed as 0; see prereg UZUPELNIENIE 1).
set -u
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s744/lock744}
O=$D/out/runs; G=$(python3 -c "import math;N=0.5+math.sqrt(1.4);print(0.1/(2*math.pi*N**3))")
for TAU in 0.5 0.25 1.0; do for SG in 1 -1; do
  for J in $(LC_ALL=C seq -f %.2f 0.30 0.02 0.70); do JR0=$J L0=1.0 GABS=$G $B 100 10 $TAU $SG 0 0 > $O/scanJ_tau${TAU}_s${SG}_Jr$J.txt; done
  for L in $(LC_ALL=C seq -f %.2f 0.80 0.02 1.20); do JR0=0.5 L0=$L GABS=$G $B 100 10 $TAU $SG 0 0 > $O/scanL_tau${TAU}_s${SG}_L$L.txt; done
done; done
echo done

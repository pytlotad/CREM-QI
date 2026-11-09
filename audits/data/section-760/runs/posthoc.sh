#!/bin/bash
# Audit 760 POST HOC extension (after ana760: zeros near 0.5 could not repeat inside the scan 0.20-1.40): J_r 1.425-1.700, kappa T 3 and soliton.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); S=${S:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s760}; O=$D/out/runs_posthoc; mkdir -p $O
G=$(python3 -c "import math;N=0.5+math.sqrt(1.4);print(0.1/(2*math.pi*N**3))")
for J in $(seq -f %.3f 1.425 0.025 1.700); do
  JR0=$J L0=1.0 GABS=$G $S/lock757 200 3 0.5 1 0 0 > $O/k3_Jr$J.txt
  SOL=1 NUT=100 KAT=0.01 GAT=0.1 JR0=$J L0=1.0 GABS=$G $S/lock759 200 0 0.5 1 0 0 > $O/sol_Jr$J.txt
done
echo done

#!/bin/bash
# Audit 746 runs (sequential), program ../section-745/lock745.cpp unchanged.  LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s746/lock745}
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/../section-745/lock745.cpp || exit 1
for E in 1e-3 -1e-3; do
  for J in $(seq -f %.2f 0.20 0.10 0.80); do for L in $(seq -f %.2f 0.70 0.10 1.30); do
    FORCE=I JR0=$J L0=$L GT=0.1 $B 300 0 0.5 1 $E 0 > $O/A_eps${E}_Jr${J}_L${L}.txt; done; done; done
for S in "0.50 2.00" "0.40 2.00" "0.60 2.00" "0.50 1.90" "0.50 2.10"; do set -- $S
  FORCE=I JR0=$1 L0=$2 GT=0.1 $B 300 0 0.5 1 1e-3 0 > $O/B_eps1e-3_Jr$1_L$2.txt; done
echo done

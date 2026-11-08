#!/bin/bash
# Audit 747 runs (sequential), program ../section-745/lock745.cpp unchanged.  LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s747/lock745}
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/../section-745/lock745.cpp || exit 1
for J in $(seq -f %.2f 0.20 0.10 0.80); do for L in $(seq -f %.2f 1.70 0.10 2.30); do
  FORCE=I JR0=$J L0=$L GT=0.1 $B 300 0 0.5 1 3.945e-4 0 > $O/C_Jr${J}_L${L}.txt; done; done
for E in 1.972e-4 7.889e-4; do for J in 0.40 0.50 0.60; do for L in 1.90 2.00 2.10; do
  FORCE=I JR0=$J L0=$L GT=0.1 $B 300 0 0.5 1 $E 0 > $O/D_eps${E}_Jr${J}_L${L}.txt; done; done; done
echo done

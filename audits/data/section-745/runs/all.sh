#!/bin/bash
# Audit 745 runs (sequential).  g++ -O3 -march=native -std=c++17 lock745.cpp.  Output -> out/runs/.  LC_ALL=C for numbers.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s745/lock745}
O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/lock745.cpp || exit 1
# control: without FORCE = lock744 (columns 1-7 vs ../section-744/out/runs/br_Jr0.5_L1.0_chi0.txt)
JR0=0.5 L0=1.0 GT=0.1 $B 300 10 0.5 1 3e-5 0 > $O/ctrl744_br_Jr0.5_L1.0_chi0.txt
for S in "0.5 1.0" "0.6 1.1"; do set -- $S
  for E in 1e-3 -1e-3 1e-4 -1e-4; do FORCE=I JR0=$1 L0=$2 GT=0.1 $B 300 10 0.5 1 $E 0 > $O/I_k10_Jr$1_L$2_eps$E.txt; done
  for E in 1e-3 -1e-3; do FORCE=I JR0=$1 L0=$2 GT=0.1 $B 300 0 0.5 1 $E 0 > $O/I_k0_Jr$1_L$2_eps$E.txt; done
done
echo done

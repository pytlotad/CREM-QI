#!/bin/bash
# Audit 748 (2): source-phase slip at EBK tori, ../section-744/lock744.cpp unchanged.  LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s748/lock744}
mkdir -p $(dirname $B); O=$D/out/slip; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/../section-744/lock744.cpp || exit 1
for T in "0.5 1.0" "0.5 2.0" "1.5 0.0" "1.5 1.0"; do set -- $T
  JR0=$1 L0=$2 GT=0.1 $B 100 10 0.5 1 0 0 > $O/slip_k10_Jr$1_L$2.txt; done
for K in 30 100; do JR0=0.5 L0=1.0 GT=0.1 $B 100 $K 0.5 1 0 0 > $O/slip_k${K}_Jr0.5_L1.0.txt; done
echo done

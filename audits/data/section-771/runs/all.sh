#!/bin/bash
# Audit 771 runs (sequential).  LC_ALL=C.  Memory cap (766).  Model: ../section-770/guide770.cpp unchanged.
set -u; export LC_ALL=C; ulimit -v 4000000
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s771/guide770}
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/../section-770/guide770.cpp || { echo BLAD kompilacja; exit 1; }
for g in 0.01:1000:50000 0.001:3000:5000000; do IFS=: read GT PER PB <<< "$g"
  for arm in m2:-1.5707963267948966 z:0; do IFS=: read a CA <<< "$arm"
    for p in A B; do for P0 in 500 $PB; do for f in 0.0 0.1 0.3 0.5 0.7 0.9 1.0; do
      GT=$GT CAUST=$CA A=1 $B $PER $f $P0 $p > $O/G_g${GT}_${a}_${p}_P${P0}_f$f.txt
    done; done; done; echo "GT $GT ramie $a"; done; done
echo ALLDONE

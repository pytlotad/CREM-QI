#!/bin/bash
# Audit 770 runs (sequential).  LC_ALL=C.  Memory cap (766).
set -u; export LC_ALL=C; ulimit -v 4000000
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s770/guide770}
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/guide770.cpp || { echo BLAD kompilacja; exit 1; }
for arm in m2:-1.5707963267948966:1 z:0:1 a0:-1.5707963267948966:0; do IFS=: read a CA AA <<< "$arm"
  for p in A B; do for P0 in 50 500; do for f in $(seq -f %.1f 0 0.1 1); do
    CAUST=$CA A=$AA $B 300 $f $P0 $p > $O/G_${a}_${p}_P${P0}_f$f.txt
  done; done; done; echo "ramie $a"; done
echo ALLDONE

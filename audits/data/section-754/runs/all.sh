#!/bin/bash
# Audit 754 run 2 (prereg754, uzupelnienie 2): bohm754.cpp, adaptive step; single-threaded processes, 4 chunks of one ensemble at a time.  LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); O=$D/out/runs; mkdir -p $O; B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s754/bohm754}
mkdir -p $(dirname $B); g++ -O3 -march=native -std=c++17 -o $B $D/bohm754.cpp || exit 1
ens(){ name=$1; N=$2; M=$3; st=$4; f=$5; c=$(( (M + 3) / 4 ))
  for i in 0 1 2 3; do m0=$((i * c)); n=$(( M - m0 < c ? M - m0 : c )); $B $N $n $st 0.25 $f $m0 > $O/${name}_p$i.txt & done; wait; }
$B 800 1 g 0.25 0.5 > $O/g_N800.txt
ens e_N800 800 2000 e 0.5
ens e_N800_f025 800 300 e 0.25
ens e_N400 400 500 e 0.5
ens e_N1600 1600 500 e 0.5
echo done

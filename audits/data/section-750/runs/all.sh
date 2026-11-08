#!/bin/bash
# Audit 750 runs (sequential).  LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s750/lock750}
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/lock750.cpp || exit 1
G=$(python3 -c "import math;N=0.5+math.sqrt(1.4);print(0.1/(2*math.pi*N**3))"); TA=$(python3 -c "import math;N=0.5+math.sqrt(1.4);print(math.pi*N**3)")
export JR0=1.5 L0=1.0 GABS=$G TAUABS=$TA
# control: without new variables lock750 = lock745 (first 7 columns of 745 file I_k0_Jr0.5_L1.0_eps1e-3)
env -u GABS -u TAUABS -u JR0 -u L0 FORCE=I JR0=0.5 L0=1.0 GT=0.1 $B 300 0 0.5 1 1e-3 0 > $O/ctrl745.txt
DRAG=1e-4 $B 1500 10 0.5 1 0 0 > $O/R0.txt
DRAG=3e-4 $B 1500 10 0.5 1 0 0 > $O/R0p.txt
DRAG=1e-4 FORCE=I $B 1500 10 0.5 1 1e-4 0 > $O/R1.txt
DRAG=1e-4 FORCE=I $B 1500 10 0.5 1 1e-3 0 > $O/R2.txt
DRAG=3e-4 FORCE=I $B 1500 10 0.5 1 1e-3 0 > $O/R2p.txt
DRAG=1e-4 FORCE=I $B 1500 0 0.5 1 1e-3 0 > $O/R3.txt
echo done

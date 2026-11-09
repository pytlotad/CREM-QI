#!/bin/bash
# Audit 769 runs (sequential).  LC_ALL=C.  Memory cap (766).
set -u; export LC_ALL=C; ulimit -v 4000000
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s769/lock769}
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/lock769.cpp || { echo BLAD kompilacja; exit 1; }
G=$(python3 -c "import math;N=0.5+math.sqrt(1.4);print(0.1/(2*math.pi*N**3))")
JR0=0.50 L0=1.0 GABS=$G $B 100 0 0.5 1 0 0 > $O/ctrl757_Jr0.50.txt
declare -A M=( [m2]=-1.5707963267948966 [z]=0 [p2]=1.5707963267948966 [m4]=-0.7853981633974483 )
for a in m2 z p2 m4; do
  for J in $(seq -f %.2f 0.20 0.05 1.70); do CAUST=${M[$a]} JR0=$J L0=1.0 GABS=$G $B 100 0 0.5 1 0 0 > $O/J_${a}_Jr$J.txt; done
  for L in $(seq -f %.2f 0.60 0.05 2.40); do CAUST=${M[$a]} JR0=0.5 L0=$L GABS=$G $B 100 0 0.5 1 0 0 > $O/L_${a}_L$L.txt; done
done
echo SKANY
for a in m2 z; do for chi in 0 1.5707963267948966; do for s in + -; do for e in 1e-5:800 1e-4:300; do E=${e%:*}; P=${e#*:}
  for S in "0.35 1.00" "0.65 1.00" "0.50 0.85" "0.50 1.15" "0.40 0.90" "0.60 1.10" "0.40 1.10" "0.60 0.90" "0.50 1.00"; do set -- $S
    CAUST=${M[$a]} JR0=$1 L0=$2 GABS=$G $B $P 0 0.5 1 ${s}$E $chi > $O/Z_${a}_c${chi:0:1}_s${s}_e${E}_Jr$1_L$2.txt
  done; done; done; done; done
echo ALLDONE

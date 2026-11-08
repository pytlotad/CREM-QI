#!/bin/bash
# Audit 744 runs (sequential).  g++ -O3 -march=native -std=c++17 lock744.cpp.  Output -> out/runs/.
set -u
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s744/lock744}
O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/lock744.cpp || exit 1
G=$(python3 -c "import math;N=0.5+math.sqrt(1.4);print(0.1/(2*math.pi*N**3))")
# H744a: kappa = 0 (gauge phase, as 743), 5 tori, 60 periods
for J in 0.30 0.40 0.50 0.60 0.70; do JR0=$J L0=1.0 GABS=$G $B 60 0 0.5 1 0 0 > $O/k0_Jr$J.txt; done
# H744b-d: slip at EBK and scans, kappa T = 10, harmony / anti-harmony, tau = T/2, T/4, T
for TAU in 0.5 0.25 1.0; do for SG in 1 -1; do
  for J in $(LC_ALL=C seq -f %.2f 0.30 0.02 0.70); do JR0=$J L0=1.0 GABS=$G $B 100 10 $TAU $SG 0 0 > $O/scanJ_tau${TAU}_s${SG}_Jr$J.txt; done
  for L in $(LC_ALL=C seq -f %.2f 0.80 0.02 1.20); do JR0=0.5 L0=$L GABS=$G $B 100 10 $TAU $SG 0 0 > $O/scanL_tau${TAU}_s${SG}_L$L.txt; done
done; done
# H744e: back-reaction
for S in "0.5 1.0" "0.6 1.1"; do set -- $S; for C in 0 1.5707963 3.1415927 4.712389; do
  JR0=$1 L0=$2 GT=0.1 $B 300 10 0.5 1 3e-5 $C > $O/br_Jr$1_L$2_chi$C.txt; done; done
echo done

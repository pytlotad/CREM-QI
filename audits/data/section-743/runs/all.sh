#!/bin/bash
# Audit 743 PDE runs (sequential).  Binary: g++ -O3 -march=native -std=c++17 field743.cpp.  Output -> out/pde/.
set -u
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s743/field743}
O=$D/out/pde; mkdir -p $O
g++ -O3 -march=native -std=c++17 -o $B $D/field743.cpp || exit 1
# H743f: 446 reproduction (default path = field446)
$B 2.3 0.6 100 0 0 0 1 3 > $O/ctrl446_eps0.txt
$B 2.3 0.6 300 3e-6 0 0 1 3 > $O/ctrl446_eps3e-6_chi0.txt
# H743c: one-period source, Gamma = 0, eps = 0
for T in "0.5 1.0 0.2" "0.6 1.1 0.2" "0.4 0.9 0.2" "0.5 1.0 0.0"; do set -- $T
  JR0=$1 L0=$2 BETA=$3 GT=0 ONEP=1 $B 5 0 0 1 > $O/onep_Jr$1_L$2_b$3.txt; done
# H743d: steady response scan, Gamma = 0.1 / T_r(0.5, 1), 40 periods
G=$(python3 -c "import math;N=0.5+math.sqrt(1.4);print(0.1/(2*math.pi*N**3))")
for J in 0.20 0.275 0.35 0.425 0.50 0.575 0.65 0.725 0.80; do
  JR0=$J L0=1.0 BETA=0.2 GABS=$G $B 40 0 0 1 > $O/scan_Jr$J.txt; done
# H743f: eps = 0 control 300 periods (GT 0.1)
JR0=0.6 L0=1.1 BETA=0.2 GT=0.1 $B 300 0 0 1 > $O/ctrl_eps0_Jr0.6_L1.1.txt
# H743e: back-reaction, gauge phase, eps 3e-5, 4 chi x 2 starts, 300 periods
for S in "0.5 1.0" "0.6 1.1"; do set -- $S
  for C in 0 1.5707963 3.1415927 4.712389; do
    JR0=$1 L0=$2 BETA=0.2 GT=0.1 $B 300 3e-5 $C 1 > $O/br_Jr$1_L$2_chi$C.txt; done; done
echo done

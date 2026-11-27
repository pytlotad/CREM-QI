#!/bin/bash
# Audit 869 (quasi-3D prototype, cause of the E1 failure): guide869.cpp; arms C (= 868 Q_01 first 300 records), S / SN (THSMOOTH, beta 0.2 / -0.05),
# D (SIGDOM, beta -0.05), T (TR2D, beta 0.2), P3 / P23 (CAUSTT -pi/3 / -2pi/3), Z0 / ZS (l = 0 starts, beta -0.05, counter / THSMOOTH).
# 3 processes, ulimit -v 2500000 each, LC_ALL=C.
set -u; export LC_ALL=C; ulimit -v 4000000
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s869; B=$S/guide869
mkdir -p $S; O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/guide869.cpp || { echo BLAD kompilacja; exit 1; }
python3 - > $D/starts869.txt <<'PY'
import random, math
random.seed(869)
for i in range(1, 11):
    ms = [(random.choice([0, 1]), 0, 0.8)]
    while len(ms) < 3:
        q = (random.randint(0, 2), random.choice([-2, -1, 1, 2]), 0.1)
        if all((q[0], q[1]) != (a, b) for a, b, _ in ms): ms.append(q)
    print(f"{i:02d} " + ';'.join(f"{a}:{b}:{w:.6f}:{random.uniform(0, 2 * math.pi):.6f}" for a, b, w in ms))
PY
echo "build $(date +%T)"; export B O
ST=$D/../section-773/starts773.txt; Z=$D/starts869.txt; P2=-1.5707963267948966; P3=-1.0471975511965976; P23=-2.0943951023931953
{ echo "C 0.2 $P2 300 01 X=0 $(awk '$1=="01"{print $2}' $ST)"
  awk -v p=$P2 '$1<=20{print "S 0.2 "p" 12000 "$1" THSMOOTH=1 "$2}' $ST
  awk -v p=$P2 '$1<=20{print "SN -0.05 "p" 12000 "$1" THSMOOTH=1 "$2}' $ST
  awk -v p=$P2 '$1<=20{print "D -0.05 "p" 12000 "$1" SIGDOM=1 "$2}' $ST
  awk -v p=$P2 '$1<=20{print "T 0.2 "p" 12000 "$1" TR2D=1 "$2}' $ST
  awk -v p=$P3 '$1<=10{print "P3 0.2 "p" 12000 "$1" X=0 "$2}' $ST
  awk -v p=$P23 '$1<=10{print "P23 0.2 "p" 12000 "$1" X=0 "$2}' $ST
  awk -v p=$P2 '{print "Z0 -0.05 "p" 12000 "$1" X=0 "$2}' $Z
  awk -v p=$P2 '{print "ZS -0.05 "p" 12000 "$1" THSMOOTH=1 "$2}' $Z; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; env $5 Q3D=1 BETA=$1 GT=0.001 CAUST=-1.5707963267948966 CAUSTT=$2 A=1 INIT="$6" $B $3 0 500 A > $O/${0}_$4.txt && echo "done ${0}_$4 $(date +%T)" || echo "FAIL ${0}_$4 $(date +%T)"'
echo "ALLDONE $(date +%T)"

#!/bin/bash
# Audit 871: local-field prototype (local871.cpp); arms K, GH, GP, GL, GZ, RH, RP, RN, S (89 runs).  3 processes, ulimit -v 2500000 each, LC_ALL=C.
set -u; export LC_ALL=C; ulimit -v 4000000
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s871; B=$S/local871
mkdir -p $S; O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/local871.cpp || { echo BLAD kompilacja; exit 1; }
echo "build $(date +%T)"; export B O
PS="M1=1 G1=2.00231930436"
{ echo "K_c KAPPA=0 N0=1.3 E0=0"; echo "K_e KAPPA=0 N0=1.0 E0=0.95"
  for n in 0.8 1.15 1.35 1.6 1.85 2.2 2.45 2.7 3.2; do for e in 0 0.5; do
    echo "GH_${n}_${e} N0=$n E0=$e"; echo "GP_${n}_${e} $PS N0=$n E0=$e"; done; echo "GL_${n} CLOCK=1 N0=$n E0=0"; done
  for n in 0.9 1.1 1.9 2.1; do echo "GZH_${n} N0=$n E0=0.95"; echo "GZP_${n} $PS N0=$n E0=0.95"; done
  for ne in 0.85,0 1.2,0 1.8,0 2.3,0 1.2,0.5 2.3,0.5; do IFS=, read n e <<< "$ne"
    echo "RH_${n}_${e} MODE=1 WDELAY=1 N0=$n E0=$e"; echo "RP_${n}_${e} MODE=1 WDELAY=1 $PS N0=$n E0=$e"; echo "RN_${n}_${e} MODE=1 WDELAY=0 N0=$n E0=$e"; done
  for ne in 1.2,0 2.3,0 1.2,0.5 2.3,0.5; do IFS=, read n e <<< "$ne"
    echo "SH0_${n}_${e} MODE=1 WDELAY=1 MS1=-0.5 MS2=0.5 N0=$n E0=$e"; echo "SH1_${n}_${e} MODE=1 WDELAY=1 MS1=0.5 MS2=0.5 N0=$n E0=$e"
    echo "SP0_${n}_${e} MODE=1 WDELAY=1 $PS MS1=-0.5 MS2=0.5 N0=$n E0=$e"; echo "SP1_${n}_${e} MODE=1 WDELAY=1 $PS MS1=0.5 MS2=0.5 N0=$n E0=$e"; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; env "$@" $B > $O/$name.txt 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

#!/bin/bash
# Audit 872: local-field prototype (../section-871/local871.cpp unchanged); starts in the capture window, TEND 6e4; arms K, GH, GP, GL, GZ,
# RH, RP, RQ, S (91 runs).  3 processes, ulimit -v 2500000 each, LC_ALL=C.
set -u; export LC_ALL=C; ulimit -v 4000000
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s872; B=$S/local871
mkdir -p $S; O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/../section-871/local871.cpp || { echo BLAD kompilacja; exit 1; }
echo "build $(date +%T)"; export B O
PS="M1=1 G1=2.00231930436"; T="TEND=60000"
{ echo "K_c KAPPA=0 N0=1.0 E0=0 $T"
  for n in 1 2 3; do for d in -0.05 -0.02 0.02 0.05; do N=$(printf '%.2f' $(echo "$n + $d" | bc -l)); for e in 0 0.5; do
    echo "GH_${N}_${e} N0=$N E0=$e $T"; echo "GP_${N}_${e} $PS N0=$N E0=$e $T"; done; done; done
  for n in 0.69 0.98 1.02 1.30 2.98 3.02; do echo "GL_${n} CLOCK=1 N0=$n E0=0 $T"; done
  for n in 0.98 1.02 1.98 2.02; do echo "GZH_${n} N0=$n E0=0.95 $T"; echo "GZP_${n} $PS N0=$n E0=0.95 $T"; done
  for ne in 0.97,0 1.03,0 1.97,0 2.03,0 1.03,0.5 2.03,0.5; do IFS=, read n e <<< "$ne"; echo "RH_${n}_${e} MODE=1 WDELAY=1 C=10 N0=$n E0=$e $T"; done
  for ne in 1.03,0 2.03,0 1.03,0.5 2.03,0.5; do IFS=, read n e <<< "$ne"
    echo "RP_${n}_${e} MODE=1 WDELAY=1 C=30 $PS N0=$n E0=$e $T"
    echo "SH0_${n}_${e} MODE=1 WDELAY=1 C=10 MS1=-0.5 MS2=0.5 N0=$n E0=$e $T"; echo "SH1_${n}_${e} MODE=1 WDELAY=1 C=10 MS1=0.5 MS2=0.5 N0=$n E0=$e $T"
    echo "SP0_${n}_${e} MODE=1 WDELAY=1 C=30 $PS MS1=-0.5 MS2=0.5 N0=$n E0=$e $T"; echo "SP1_${n}_${e} MODE=1 WDELAY=1 C=30 $PS MS1=0.5 MS2=0.5 N0=$n E0=$e $T"; done
  for n in 1.03 2.03; do echo "RQ_${n}_0 MODE=1 WDELAY=1 C=10 $PS N0=$n E0=0 $T"; done; } > $S/jobs.txt
echo "jobs $(wc -l < $S/jobs.txt)"
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; env "$@" $B > $O/$name.txt 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"' < $S/jobs.txt
echo "ALLDONE $(date +%T)"

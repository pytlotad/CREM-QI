#!/bin/bash
# Audit 873: n = 1 l ~ 0 fall -- kernel width w / grid dx / eccentricity scan (../section-871/local871.cpp unchanged); 16 runs.
set -u; export LC_ALL=C; ulimit -v 4000000
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s873; B=$S/local871
mkdir -p $S; O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/../section-871/local871.cpp || { echo BLAD kompilacja; exit 1; }
echo "build $(date +%T)"; export B O
PS="M1=1 G1=2.00231930436"; T="TEND=60000"
{ echo "K KAPPA=0 N0=1.02 E0=0.95 $T"
  for n in 0.98 1.02; do echo "D_H_${n} N0=$n E0=0.95 W=0.25 DX=0.025 $T"; done
  for wd in 0.05,0.01,W05 0.10,0.02,W10; do IFS=, read w d tag <<< "$wd"; for n in 0.98 1.02; do
    echo "${tag}_H_${n} N0=$n E0=0.95 W=$w DX=$d $T"; echo "${tag}_P_${n} $PS N0=$n E0=0.95 W=$w DX=$d $T"; done; done
  for e in 0.80 0.85 0.90; do echo "E_H_${e} N0=1.02 E0=$e $T"; done
  echo "X_H_2.02 N0=2.02 E0=0.95 W=0.5 DX=0.1 $T"; echo "X_P_2.02 $PS N0=2.02 E0=0.95 W=1.0 DX=0.2 $T"; } > $S/jobs.txt
echo "jobs $(wc -l < $S/jobs.txt)"
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; env "$@" $B > $O/$name.txt 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"' < $S/jobs.txt
echo "ALLDONE $(date +%T)"

#!/bin/bash
# Audit 875: local field + Landau-Lifshitz + LW at c = 10, 30, 137; 32 runs + byte control (local875 LL=0 vs local871).
set -u; export LC_ALL=C; ulimit -v 4000000
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s875; B=$S/local875
mkdir -p $S; O=$D/out/runs; mkdir -p $O
g++ -O3 -march=native -std=c++17 -o $B $D/local875.cpp || { echo BLAD kompilacja; exit 1; }
g++ -O3 -march=native -std=c++17 -o $S/local871 $D/../section-871/local871.cpp || { echo BLAD kompilacja 871; exit 1; }
echo "build $(date +%T)"; export B O
# byte control (H875c)
( ulimit -v 2500000
  env MODE=1 WDELAY=1 C=10 N0=1.02 E0=0 TEND=3000 $S/local871 > $O/C1_871.txt; env MODE=1 WDELAY=1 C=10 N0=1.02 E0=0 TEND=3000 $B > $O/C1_875.txt
  env N0=1.02 E0=0.95 TEND=3000 $S/local871 > $O/C2_871.txt; env N0=1.02 E0=0.95 TEND=3000 $B > $O/C2_875.txt )
for k in 1 2; do cmp -s $O/C${k}_871.txt $O/C${k}_875.txt && echo "C$k IDENT" || echo "C$k ROZNE"; done
echo "control $(date +%T)"
python3 $D/runs/jobs875.py > $S/jobs.txt; echo "jobs $(wc -l < $S/jobs.txt)"; cat $S/jobs.txt
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; env "$@" $B > $O/$name.txt 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"' < $S/jobs.txt
echo "ALLDONE $(date +%T)"

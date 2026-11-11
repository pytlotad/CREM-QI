#!/bin/bash
# Audit 775 runs: beta 0.05 and 0.8, starts 01-20 of ../section-773/starts773.txt; 3 parallel (xargs -P 3), ulimit -v 4 GB each.  LC_ALL=C.  Model ../section-772/guide772.cpp unchanged.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s775/guide772}
mkdir -p $(dirname $B); O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/../section-772/guide772.cpp || { echo BLAD kompilacja; exit 1; }
export B O
( ulimit -v 4000000; BETA=0.2 GT=0.0003 CAUST=-1.5707963267948966 A=1 INIT="$(awk '$1=="01"{print $2}' $D/../section-773/starts773.txt)" $B 300 0 500 A > $O/C_ctrl.txt ) && echo "done C"
{ for a in "B05 0.05" "B80 0.8"; do awk -v a="$a" '$1<=20{print a" "$1" "$2}' $D/../section-773/starts773.txt; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 4000000; BETA=$1 GT=0.0003 CAUST=-1.5707963267948966 A=1 INIT="$3" $B 20000 0 500 A > $O/${0}_$2.txt && echo "done ${0}_$2"'
echo ALLDONE

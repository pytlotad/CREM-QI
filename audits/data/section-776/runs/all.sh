#!/bin/bash
# Audit 776: drive spectra on frozen start tori; 3 parallel (xargs -P 3), ulimit -v 4 GB each.  LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s776/predict776}
mkdir -p $(dirname $B); O=$D/out/pred; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/predict776.cpp || { echo BLAD kompilacja; exit 1; }
export B O
{ for c in "c1 -1.5707963267948966" "c0 0"; do
    awk -v c="$c" '{print "S772_b02 0.2 "c" "$1" "$2}' $D/../section-772/starts772.txt
    awk -v c="$c" '$1<=40{print "S773_b02 0.2 "c" "$1" "$2}' $D/../section-773/starts773.txt
    awk -v c="$c" '$1<=20{print "S773_b005 0.05 "c" "$1" "$2}' $D/../section-773/starts773.txt
    awk -v c="$c" '$1<=20{print "S773_b08 0.8 "c" "$1" "$2}' $D/../section-773/starts773.txt; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 4000000; BETA=$1 CAUST=$3 INIT="$5" $B 300 > $O/${0}_${2}_$4.txt && echo "done ${0}_${2}_$4"'
echo ALLDONE

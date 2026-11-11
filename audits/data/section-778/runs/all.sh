#!/bin/bash
# Audit 778: (1) code control, (2) reruns to k_d + 20 with DUMP, (3) drive spectra on the current torus at k_d-50, k_d-10, k_d.  3 parallel, ulimit -v 4 GB each.  LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); SP=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s778
mkdir -p $SP; O=$D/out/runs; Q=$D/out/pred; mkdir -p $O $Q
g++ -O3 -march=native -std=c++17 -o $SP/guide778 $D/guide778.cpp && g++ -O3 -march=native -std=c++17 -o $SP/predict778 $D/predict778.cpp || { echo BLAD kompilacja; exit 1; }
export SP O Q
( ulimit -v 4000000; BETA=0.2 GT=0.0003 CAUST=-1.5707963267948966 A=1 INIT="$(awk '$1=="01"{print $2}' $D/../section-773/starts773.txt)" $SP/guide778 300 0 500 A > $O/C_ctrl.txt ) && echo "done C"
xargs -P 3 -L 1 bash -c 'ulimit -v 4000000; BETA=$2 GT=0.0003 CAUST=-1.5707963267948966 A=1 INIT="$4" DUMP=$O/${0}_$1.dump $SP/guide778 $(($3 + 20)) 0 500 A > $O/${0}_$1.txt && echo "done ${0}_$1"' < $D/jobs778.txt
while read a i b kd ini; do for off in 50 10 0; do k=$((kd - off)); [ $k -ge 1 ] || continue
  read J L < <(awk -v k=$k '!/^#/ && $1==k {print $3, $4}' $O/${a}_$i.txt); echo "$a $i $b $off $J $L"; done; done < $D/jobs778.txt |
xargs -P 3 -L 1 bash -c 'ulimit -v 4000000; BETA=$2 CAUST=-1.5707963267948966 JT=$4 LT=$5 INIT="0:0:1:0" $SP/predict778 300 > $Q/${0}_$1_m$3.txt && echo "done P ${0}_$1_m$3"'
echo ALLDONE

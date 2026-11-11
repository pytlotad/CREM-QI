#!/bin/bash
# Audit 779: one-period drive integrals I_n on frozen tori -- (1) tori of the 388 existing drive spectra (776, 778), (2) EBK tori n_r 0..8, |m| <= 5 at beta 0.05/0.2/0.8, CAUST -pi/2 and 0.
# No model runs.  LC_ALL=C, ulimit -v 4 GB.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); SP=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s779
mkdir -p $SP $D/out
g++ -O3 -march=native -std=c++17 -o $SP/overlap779 $D/overlap779.cpp || { echo BLAD kompilacja; exit 1; }
for s in 776 778; do for f in $D/../section-$s/out/pred/*.txt; do
  head -1 $f | awk -v k="$s:$(basename $f .txt)" '{for(i=1;i<=NF;i++){if($i=="beta")b=$(i+1); if($i=="caust")c=$(i+1); if($i=="J0")j=$(i+1); if($i=="L0")l=$(i+1)}; print k, b, c, j, l}'
done; done > $D/out/tori_spectra.txt
for b in 0.05 0.2 0.8; do for c in -1.5707963267948966 0; do for nr in $(seq 0 8); do for m in $(seq -5 5); do
  echo "ebk:$b:$c:$nr:$m $b $c $(awk -v n=$nr 'BEGIN{printf "%.1f", n + 0.5}') $m"; done; done; done; done > $D/out/tori_ebk.txt
( ulimit -v 4000000; $SP/overlap779 < $D/out/tori_spectra.txt > $D/out/I_spectra.txt ) && echo "done spectra"
( ulimit -v 4000000; $SP/overlap779 < $D/out/tori_ebk.txt > $D/out/I_ebk.txt ) && echo "done ebk"
echo ALLDONE

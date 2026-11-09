#!/bin/bash
# Audit 767: 28 states, one process per state, sequential, memory-capped; states already done are skipped.
export LC_ALL=C OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=/home/teddy/Projekty/CREM-QI/audits/data/section-767; cd $D
for n in 2 3 4 5 6 7 8; do for ((l=1; l<n; l++)); do
  [ -f out/stany/${n}_${l}.txt ] && continue
  ( ulimit -v 4000000; nice -n 10 python3 stan767.py $n $l ) || echo "BLAD $n $l (kod $?)"
done; done
echo ALLDONE

#!/bin/bash
S="$1"; cd /home/teddy/Projekty/Positronium
export ONLY_PARA=1 CREM_L_BALANCE=1 CREM_EMISSION_REACH=1 CREM_MICROCANONICAL_START=1
for w in 0 1 2 3; do
  ( for i in $(seq $w 4 23); do $S/sec329/cause $i > $S/sec329/traj-$i.txt 2>&1; done ) &
done
wait
touch $S/sec329/DONE

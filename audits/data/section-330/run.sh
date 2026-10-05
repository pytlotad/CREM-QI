#!/bin/bash
S="$1"; cd /home/teddy/Projekty/Positronium
make validation-small > "$S/sec330/walidacja-small.txt" 2>&1
export ONLY_PARA=1 CREM_L_BALANCE=1 CREM_EMISSION_REACH=1 CREM_MICROCANONICAL_START=1
for w in 0 1 2 3; do ( for i in $(seq $w 4 23); do "$S/sec330/cause" $i > "$S/sec330/traj-$i.txt" 2>&1; done ) & done
wait
unset ONLY_PARA CREM_L_BALANCE CREM_EMISSION_REACH CREM_MICROCANONICAL_START
g++ -std=c++20 -O3 -march=native -I . $(root-config --cflags) tools/para_ortho_lifetimes.cpp -o "$S/sec330/pol" $(root-config --libs) > "$S/sec330/pol-build.txt" 2>&1
"$S/sec330/pol" 48 42 1200 det 1 nofloor branch > "$S/sec330/A-kolo.txt" 2>&1
CREM_MICROCANONICAL_START=1 "$S/sec330/pol" 48 42 1200 det 1 nofloor branch > "$S/sec330/B-zespol.txt" 2>&1
git checkout -- distributions/ 2>/dev/null
touch "$S/sec330/DONE"

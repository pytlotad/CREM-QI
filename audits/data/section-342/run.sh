#!/bin/bash
S="$1"; cd /home/teddy/Projekty/Positronium
make validation-small > "$S/sec342/walidacja-small.txt" 2>&1
g++ -std=c++20 -O3 -march=native -I . $(root-config --cflags) tools/para_ortho_lifetimes.cpp -o "$S/sec342/pol" $(root-config --libs) > "$S/sec342/pol-build.txt" 2>&1
"$S/sec342/pol" 48 42 1200 det 1 nofloor branch > "$S/sec342/A-kolo-swobodne.txt" 2>&1
CREM_SPIN_QUANTIZATION=1 "$S/sec342/pol" 48 42 1200 det 1 nofloor branch > "$S/sec342/A-kolo-kwant.txt" 2>&1
CREM_MICROCANONICAL_START=1 "$S/sec342/pol" 48 42 1200 det 1 nofloor branch > "$S/sec342/B-zespol-swobodne.txt" 2>&1
CREM_MICROCANONICAL_START=1 CREM_SPIN_QUANTIZATION=1 "$S/sec342/pol" 48 42 1200 det 1 nofloor branch > "$S/sec342/B-zespol-kwant.txt" 2>&1
git checkout -- distributions/ 2>/dev/null
touch "$S/sec342/DONE"

#!/bin/bash
S="$1"; cd /home/teddy/Projekty/Positronium
make validation-small > "$S/sec328/walidacja-small.txt" 2>&1
g++ -std=c++20 -O3 -march=native -I . $(root-config --cflags) tools/para_ortho_lifetimes.cpp -o "$S/sec328/pol" $(root-config --libs) > "$S/sec328/pol-build.txt" 2>&1
"$S/sec328/pol" 48 42 1200 det 1 nofloor branch > "$S/sec328/A-kolo.txt" 2>&1
CREM_MICROCANONICAL_START=1 "$S/sec328/pol" 48 42 1200 det 1 nofloor branch > "$S/sec328/B-zespol.txt" 2>&1
touch "$S/sec328/DONE"

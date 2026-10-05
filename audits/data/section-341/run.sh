#!/bin/bash
S="$1"; cd /home/teddy/Projekty/Positronium
make validation-small > "$S/sec341/walidacja-small.txt" 2>&1
g++ -std=c++20 -O3 -march=native -I . $(root-config --cflags) tools/para_ortho_lifetimes.cpp -o "$S/sec341/pol" $(root-config --libs) > /dev/null 2>&1
"$S/sec341/pol" 48 42 1200 det 1 nofloor branch > "$S/sec341/A-kolo.txt" 2>&1
git checkout -- distributions/ 2>/dev/null
touch "$S/sec341/DONE"

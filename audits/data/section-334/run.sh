#!/bin/bash
S="$1"; cd /home/teddy/Projekty/Positronium
make validation-small > "$S/sec334/walidacja-small.txt" 2>&1
g++ -std=c++20 -O3 -march=native -I . $(root-config --cflags) tools/para_ortho_lifetimes.cpp -o "$S/sec334/pol" $(root-config --libs) > "$S/sec334/pol-build.txt" 2>&1
"$S/sec334/pol" 48 42 1200 det 1 nofloor branch > "$S/sec334/A-kolo.txt" 2>&1
git checkout -- distributions/ 2>/dev/null
touch "$S/sec334/DONE"

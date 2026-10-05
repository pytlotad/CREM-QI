#!/bin/bash
S="$1"; cd /home/teddy/Projekty/Positronium
make validation-small > "$S/sec340/walidacja-small.txt" 2>&1
g++ -O3 -march=native -std=c++20 -I. $(root-config --cflags) "$S/sec338/cause3.cpp" -o "$S/sec340/cause" $(root-config --libs) > /dev/null 2>&1
( "$S/sec340/cause" 7 > "$S/sec340/c7.out" 2>/dev/null ) &
( "$S/sec340/cause" 20 > "$S/sec340/c20.out" 2>/dev/null ) &
wait
g++ -std=c++20 -O3 -march=native -I . $(root-config --cflags) tools/para_ortho_lifetimes.cpp -o "$S/sec340/pol" $(root-config --libs) > /dev/null 2>&1
"$S/sec340/pol" 48 42 1200 det 1 nofloor branch > "$S/sec340/A-kolo.txt" 2>&1
git checkout -- distributions/ 2>/dev/null
touch "$S/sec340/DONE"

#!/bin/bash
S="$1"; cd /home/teddy/Projekty/Positronium
while [ ! -e "$S/lev/DONE" ]; do sleep 30; done
timeout 14000 ./positronium --mode statistical --phenomenon 6 --runs 30 --seed 42 > "$S/sec320/n1-30par.txt" 2>&1
touch "$S/sec320/DONE"

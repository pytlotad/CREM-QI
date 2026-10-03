#!/bin/bash
S="$1"; cd /home/teddy/Projekty/Positronium
while pgrep -x positronium >/dev/null; do sleep 30; done
timeout 14000 ./positronium --mode statistical --phenomenon 6 --level 2 --contact-j0 0.035 --runs 30 --seed 42 > "$S/lev/n2.txt" 2>&1
timeout 30000 ./positronium --mode statistical --phenomenon 6 --level 3 --contact-j0 0.0233333 --runs 30 --seed 42 > "$S/lev/n3.txt" 2>&1
touch "$S/lev/DONE"

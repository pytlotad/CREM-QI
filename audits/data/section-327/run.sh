#!/bin/bash
S="$1"; cd /home/teddy/Projekty/Positronium
$S/sec326/pol 48 42 1200 det 1 nofloor x branch > "$S/sec327/A-kolo.txt" 2>&1
CREM_MICROCANONICAL_START=1 $S/sec326/pol 48 42 1200 det 1 nofloor x branch > "$S/sec327/B-zespol.txt" 2>&1
touch "$S/sec327/DONE"

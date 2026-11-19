#!/bin/bash
# Audit 811: probe p811 on tree d277071b, new and off (CREM_NO_SOMMERFELD_EXACT=1).  LC_ALL=C, memory cap (766).
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s811; D=$R/audits/data/section-811
rm -rf $S; mkdir -p $S/t; ( cd $R && git archive d277071b | tar -x -C $S/t )
( cd $S/t && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules $D/p811.cpp -o $S/p811 $(root-config --libs) ) || { echo BLAD; exit 1; }
$S/p811 > $D/out/p811_new.txt && CREM_NO_SOMMERFELD_EXACT=1 $S/p811 > $D/out/p811_off.txt && echo "ALLDONE $(date +%T)"

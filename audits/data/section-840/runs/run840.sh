#!/bin/bash
# Audit 840: probe p840 (hydrogen test particle vs Dirac FW classical limit) -- F0 / F1 x 5 nodes x {base, N2, h2} + control CL.
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s840; D=$R/audits/data/section-840; O=$D/out/runs
mkdir -p $S $O
( cd $R && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-840/p840.cpp -o $S/p840 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"; export S O
declare -A LF=([u5a]=0.5120091953982716 [u5b]=0.5652161499007791 [u5c]=0.66666666666666663 [u5d]=0.8125025652578265 [u5e]=0.9551918755393795)
grep -E '^H[01] ' $D/out/params839.txt | while read sk nd pair rp e a b c x y z; do
  for v in "base 4096 65536" "N2 8192 65536" "h2 4096 131072"; do set -- $v; echo "$sk $nd $1 ${LF[$nd]} $x $y $z $2 $3 R"; done
  echo "$sk $nd CL ${LF[$nd]} $x $y $z 1024 65536 C"
done | xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; sk=$0; nd=$1; var=$2; shift 2
  "$S/p840" "$@" > $O/${sk}_${nd}_${var}.txt 2>&1 && { [ "$7" = C ] && CREM_THOMAS_TULCZYJEW_INSTANT=1 "$S/p840" "$@" >> $O/${sk}_${nd}_${var}.txt 2>&1; true; } \
  && echo "done $sk $nd $var $(date +%T)" || echo "FAIL $sk $nd $var $(date +%T)"'
echo "ALLDONE $(date +%T)"

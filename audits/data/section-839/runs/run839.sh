#!/bin/bash
# Audit 839: probe p839 at the 15 (system, node) orbits of 838 -- base (N 4096, h = T/65536), N 8192, h = T/131072; control C (L and T).
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s839; D=$R/audits/data/section-839; O=$D/out/runs
mkdir -p $S $O
( cd $R && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-839/p839.cpp -o $S/p839 $(root-config --libs) ) || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"; export S O
while read sk nd pair rp e a b c d f g; do
  for v in "base 4096 65536" "N2 8192 65536" "h2 4096 131072"; do set -- $v
    echo "$sk $nd $1 $2 $3 $pair $rp $e $a $b $c $d $f $g R"; done
  echo "$sk $nd CL 4096 65536 $pair $rp $e $a $b $c $d $f $g C"
done < $D/out/params839.txt |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; sk=$0; nd=$1; var=$2; N=$3; hd=$4; shift 4
  "$S/p839" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" $N $hd "${10}" > $O/${sk}_${nd}_${var}.txt 2>&1 \
  && { [ "${10}" = C ] && CREM_THOMAS_TULCZYJEW_INSTANT=1 "$S/p839" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" $N $hd C >> $O/${sk}_${nd}_${var}.txt 2>&1; true; } \
  && echo "done $sk $nd $var $(date +%T)" || echo "FAIL $sk $nd $var $(date +%T)"'
echo "ALLDONE $(date +%T)"

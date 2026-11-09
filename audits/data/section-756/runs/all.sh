#!/bin/bash
# Audit 756: p585 probe (unchanged) on the engine tree 0f9406c4; print-only.  Sequential.  LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); R=$(cd $D/../../.. && pwd); O=$D/out/runs; mkdir -p $O
B=${B:-/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s756/p585}
( cd $R && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-585/p585.cpp -o $B $(root-config --libs) ) || exit 1
for pair in H Ps; do
  P=""; [ $pair = H ] && P="PAIR=proton,electron"
  for s in "2.15 1.5" "2.3 1.5" "2.5 1.5" "2.7 1.5" "2.85 1.5" "2 1.1" "2 1.3" "2 1.7" "2 1.9" "3.4 2.2"; do
    set -- $s; n=$1; L=$2; f=$(python3 -c "print(repr($L/$n))")
    for w in 1 0; do
      W=""; [ $w = 0 ] && W="CREM_WAVE_LOCK=0"
      ( cd $O && env CREM_EMISSION_LEDGER=1 CREM_DEBUG_LADDER=1 $P $W NLEV=$n LFRAC=$f $B 0 > ${pair}_n${n}_L${L}_W$w.out 2> ${pair}_n${n}_L${L}_W$w.err )
    done
  done
done
echo done

#!/bin/bash
S="$1"; cd "$S/sec339"
export CREM_E1_CORRECTION=1 CREM_L_BALANCE=1 CREM_EMISSION_REACH=1
for w in 0 1 2 3; do
  ( for i in $(seq $w 4 47); do
      ONLY_PARA=1  "$S/sec338/cause" $i > p$i.out 2> p$i.err
      ONLY_ORTHO=1 "$S/sec338/cause" $i > o$i.out 2> o$i.err
    done ) &
done
wait; touch DONE

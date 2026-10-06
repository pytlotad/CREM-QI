#!/bin/bash
# job.sh LEVEL J0 VARIANT IDX -> out/nLEVEL-lJ0-VARIANT-IDX.txt ; VARIANT A (reflection) or C (Delta l = +1 below hbar)
cd "$(dirname "$0")"; n=$1; j=$2; v=$3; i=$4; o=out/n$n-j$j-$v-$i; [ -f $o.txt ] && exit 0
export LEVEL=$n CREM_INITIAL_ANGULAR_MOMENTUM=$j CREM_EMISSION_REACH=1
[ $v = C ] && export CREM_DL_PLUS_BELOW_HBAR=1
./langer $i > $o.tmp && mv $o.tmp $o.txt

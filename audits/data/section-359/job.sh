#!/bin/bash
cd "$(dirname "$0")"; v=$1; i=$2; [ -f out/$v-$i.txt ] && exit 0
export CREM_EMISSION_REACH=1
[ $v = C ] && export CREM_DL_PLUS_BELOW_HBAR=1
./trace2 $i > out/$v-$i.tmp && mv out/$v-$i.tmp out/$v-$i.txt

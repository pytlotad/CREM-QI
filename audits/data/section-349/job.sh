#!/bin/bash
# job.sh VARIANT IDX -> out/VARIANT-IDX.txt ; A = action rule, C = action rule + Delta l = +1 below hbar
cd "$(dirname "$0")"; v=$1; i=$2; [ -f out/$v-$i.txt ] && exit 0
export CREM_EMISSION_REACH=1 CREM_ACTION_PHOTON=1
case $v in C) export CREM_DL_PLUS_BELOW_HBAR=1;; esac
./trace2 $i > out/$v-$i.tmp && mv out/$v-$i.tmp out/$v-$i.txt

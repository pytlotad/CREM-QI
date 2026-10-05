#!/bin/bash
# job.sh VARIANT IDX -> out/VARIANT-IDX.txt ; R = 343 set at n = 1 (production), A/B/C = n = 2 variants
cd "$(dirname "$0")"; v=$1; i=$2; [ -f out/$v-$i.txt ] && exit 0; export CREM_EMISSION_REACH=1; bin=./trace2
case $v in B) export CREM_CLOSE_BELOW_HBAR=1;; C) export CREM_DL_PLUS_BELOW_HBAR=1;; R) bin=./trace1;; esac
$bin $i > out/$v-$i.tmp && mv out/$v-$i.tmp out/$v-$i.txt

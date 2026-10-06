#!/bin/bash
cd "$(dirname "$0")"; i=$1; [ -f out/$i.txt ] && exit 0
CREM_ACTION_PHOTON=1 CREM_HOLD_AFTER_CLOSURE=5e-10 ./trace2 $i > out/$i.tmp && mv out/$i.tmp out/$i.txt

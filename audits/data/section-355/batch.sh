#!/bin/bash
cd "$(dirname "$0")"
run(){ i=$1; LEVEL=2 CREM_INITIAL_ANGULAR_MOMENTUM=0.75 CREM_EMISSION_REACH=1 ./langer $i > out/$i.tmp && mv out/$i.tmp out/$i.txt; }
export -f run 2>/dev/null
for i in 0 1 2 3 4 5; do echo $i; done | xargs -P 4 -I{} bash -c 'cd "'"$PWD"'"; LEVEL=2 CREM_INITIAL_ANGULAR_MOMENTUM=0.75 CREM_EMISSION_REACH=1 ./langer {} > out/{}.tmp && mv out/{}.tmp out/{}.txt'
LEVEL=2 CREM_INITIAL_ANGULAR_MOMENTUM=0.75 CREM_HOLD_AFTER_CLOSURE=5e-11 stdbuf -oL ./langer 0 > hold0.txt
echo BATCH-DONE

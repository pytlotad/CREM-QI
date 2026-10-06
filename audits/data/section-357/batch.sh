#!/bin/bash
cd "$(dirname "$0")"
for i in 0 1 2 3 4 5; do echo $i; done | xargs -P 3 -I{} sh -c 'LEVEL=2 CREM_INITIAL_ANGULAR_MOMENTUM=0.75 CREM_EMISSION_REACH=1 ./langer {} > out/{}.tmp && mv out/{}.tmp out/{}.txt' &
LEVEL=2 CREM_INITIAL_ANGULAR_MOMENTUM=0.75 CREM_HOLD_AFTER_CLOSURE=2e-9 stdbuf -oL ./langer_o 0 > hold2ns_o.txt 2>&1
wait; echo BATCH-DONE

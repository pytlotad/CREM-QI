#!/bin/bash
S="$1"; cd /home/teddy/Projekty/Positronium
for sm in 0.30 0.15; do (PROBE_SMAX=$sm ONLY_PARA=1 CREM_MICROCANONICAL_START=1 $S/sec328/cause-smax 14 > $S/sec328/v2-idx14-s$sm.txt 2>&1) & done
for sm in 0.30 0.15; do (PROBE_SMAX=$sm ONLY_PARA=1 $S/sec328/cause-smax 7 > $S/sec328/v2-circ7-s$sm.txt 2>&1) & done
wait
(CREM_CONTINUOUS_ORBIT_CREDIT=1 PROBE_SMAX=0.30 ONLY_PARA=1 CREM_MICROCANONICAL_START=1 $S/sec328/cause-smax 14 > $S/sec328/v2-idx14-old.txt 2>&1) &
(CREM_STOCHASTIC_SMAX_CAP=1 PROBE_SMAX=0.30 ONLY_PARA=1 $S/sec328/cause-smax 7 > $S/sec328/v2-circ7-cap.txt 2>&1) &
wait
for f in v2-idx14-s0.30 v2-idx14-s0.15 v2-idx14-old v2-circ7-s0.30 v2-circ7-s0.15 v2-circ7-cap; do echo "$f: $(cat $S/sec328/$f.txt)"; done

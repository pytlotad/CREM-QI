#!/bin/bash
S="$1"; cd /home/teddy/Projekty/Positronium
(PROBE_SMAX=0.30 ONLY_PARA=1 CREM_MICROCANONICAL_START=1 $S/sec328/cause-smax 14 > $S/sec328/idx14-s030.txt 2>&1) &
(PROBE_SMAX=0.15 ONLY_PARA=1 CREM_MICROCANONICAL_START=1 $S/sec328/cause-smax 14 > $S/sec328/idx14-s015.txt 2>&1) &
(CREM_CONTINUOUS_ORBIT_CREDIT=1 PROBE_SMAX=0.30 ONLY_PARA=1 CREM_MICROCANONICAL_START=1 $S/sec328/cause-smax 14 > $S/sec328/idx14-old.txt 2>&1) &
(timeout 1200 ./positronium --mode statistical --phenomenon 1 --level 1 --radiation-reaction stochastic --runs 1 --seed 7 --crem-wallclock-budget-s 300 > $S/sec328/circ-328.txt 2>&1) &
wait
git checkout -- distributions/ 2>/dev/null
for f in idx14-s030 idx14-s015 idx14-old; do echo "$f: $(cat $S/sec328/$f.txt)"; done
diff <(grep -v "^CREM calib" $S/sec325/circ-after.txt) <(grep -v "^CREM calib" $S/sec328/circ-328.txt) | head -30

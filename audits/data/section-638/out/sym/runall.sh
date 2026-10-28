#!/bin/bash
cd /tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r638/sym
{ for sys in ps h; do for B in 0 0.1; do for T in id P C CP; do echo "./t1.sh $sys $B $T"; done; done; done
  for cs in pps ops free h; do for B in 0 0.1; do for T in id P C CP; do echo "BUDGET=3000 ./t2.sh $cs $B $T"; done; done; done; } | xargs -P 2 -I{} bash -c "{}"
python3 /home/teddy/Projekty/CREM-QI/audits/data/section-623/compare623.py /tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r638/sym/T1 > /tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r638/sym/T1_compare.txt 2>&1
python3 /home/teddy/Projekty/CREM-QI/audits/data/section-623/compare623b.py /tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r638/sym/T2 > /tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r638/sym/T2_compare.txt 2>&1
echo DONE > /tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r638/sym/done.txt

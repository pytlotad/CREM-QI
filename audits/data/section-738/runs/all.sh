#!/bin/bash
# sequential: build p585 (base tree 3cf1b164 and new tree), probes old/new/new_off, regression + validation
S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s738; D=/home/teddy/Projekty/CREM-QI/audits/data/section-738/runs
mkdir -p $S/base; ( cd /home/teddy/Projekty/CREM-QI && git archive 3cf1b164 | tar -x -C $S/base )
( cd $S/base && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-585/p585.cpp -o $S/p585_base $(root-config --libs) )
cd /home/teddy/Projekty/CREM-QI && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-585/p585.cpp -o $S/p585_new $(root-config --libs)
bash $D/probes.sh $S/p585_base $S/probes/old
bash $D/probes.sh $S/p585_new $S/probes/new
bash $D/probes.sh $S/p585_new $S/probes/new_off CREM_EBK_ORDER0=1
for t in H2p H3p H4p H3s H3d Ps2p; do s=SAME; cmp -s $S/probes/old/$t.out $S/probes/new_off/$t.out && cmp -s $S/probes/old/$t.err $S/probes/new_off/$t.err || s=DIFF; echo "probe $t old vs new_off $s"; done > $S/probes_cmp.txt
for t in H2p H3p H4p H3s H3d Ps2p; do s=SAME; cmp -s $S/probes/old/$t.err $S/probes/new/$t.err || s=DIFF; echo "probe $t old vs new $s"; done >> $S/probes_cmp.txt
bash $D/reg.sh
echo ALLDONE

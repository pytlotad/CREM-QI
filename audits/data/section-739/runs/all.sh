#!/bin/bash
# sequential: build p585 and p623 (base tree 637edfd9 and new tree), probes old/new/new_off, H739d statistics (p623 i = 0..63, old/new),
# regression + validation
S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s739; D=/home/teddy/Projekty/CREM-QI/audits/data/section-739/runs
mkdir -p $S/base; ( cd /home/teddy/Projekty/CREM-QI && git archive 637edfd9 | tar -x -C $S/base )
for p in 585 623; do
  ( cd $S/base && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-$p/p$p.cpp -o $S/p${p}_base $(root-config --libs) )
  ( cd /home/teddy/Projekty/CREM-QI && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-$p/p$p.cpp -o $S/p${p}_new $(root-config --libs) )
done
bash $D/probes.sh $S/p585_base $S/probes/old
bash $D/probes.sh $S/p585_new $S/probes/new
bash $D/probes.sh $S/p585_new $S/probes/new_off CREM_CLASSICAL_S_TRIPLET=1
for t in H2p H3p H4p H3s H3d Ps2p; do s=SAME; cmp -s $S/probes/old/$t.out $S/probes/new_off/$t.out && cmp -s $S/probes/old/$t.err $S/probes/new_off/$t.err || s=DIFF; echo "probe $t old vs new_off $s"; done > $S/probes_cmp.txt
for t in H2p H3p H4p H3s H3d Ps2p; do s=SAME; cmp -s $S/probes/old/$t.out $S/probes/new/$t.out && cmp -s $S/probes/old/$t.err $S/probes/new/$t.err || s=DIFF; echo "probe $t old vs new $s"; done >> $S/probes_cmp.txt
mkdir -p $S/stat
for v in base new; do
  ( cd $S/stat && env MODE=est PAIR=proton,electron PHEN=1 QSPIN=1 NLEV=2 LFRAC=0.75 BUDGET=3000 CREM_EMISSION_LEDGER=1 \
      $S/p623_$v $(seq 0 63) > $v.out 2> $v.err )
done
bash $D/reg.sh
echo ALLDONE

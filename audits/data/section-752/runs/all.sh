#!/bin/bash
# Audit 752: byte-identity of the descriptive change (modeBeatEnergy).  Base tree 5ed7ea44 (751) vs working tree.  Sequential.
set -u; export LC_ALL=C
S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s752; D=/home/teddy/Projekty/CREM-QI/audits/data/section-752
rm -rf $S; mkdir -p $S/old $S/new
( cd /home/teddy/Projekty/CREM-QI && git archive 5ed7ea44 | tar -x -C $S/old ); ( cd /home/teddy/Projekty/CREM-QI && rsync -a --exclude .git --exclude audits ./ $S/new/ )
mkdir -p $S/new/audits/data/section-585; cp /home/teddy/Projekty/CREM-QI/audits/data/section-585/p585.cpp $S/new/audits/data/section-585/
for v in old new; do ( cd $S/$v && g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-585/p585.cpp -o $S/p585_$v $(root-config --libs) ); done
for v in old new; do bash /home/teddy/Projekty/CREM-QI/audits/data/section-741/runs/probes.sh $S/p585_$v $S/probes/$v; done
for t in H2p H3p H4p H3s H3d Ps2p; do s=SAME; cmp -s $S/probes/old/$t.out $S/probes/new/$t.out && cmp -s $S/probes/old/$t.err $S/probes/new/$t.err || s=DIFF; echo "probe $t old vs new $s"; done > $D/out/probes752.txt
( cd $S/old && make -j4 crem_interactions >/dev/null 2>&1 ); ( cd $S/new && make -j4 crem_interactions crem_interactions_validation >/dev/null 2>&1 )
( cd $S/new && ./crem_interactions_validation > $D/out/validation.txt 2>&1; echo EXIT $? >> $D/out/validation.txt )
run(){ v=$1; tag=$2; shift 2; mkdir -p $S/reg/$v/$tag; ( cd $S/reg/$v/$tag && $S/$v/crem_interactions "$@" > out.txt 2> err.txt ); }
for v in old new; do
run $v 1 --mode statistical --phenomenon 1 --runs 16 --seed 466
run $v 2 --mode statistical --phenomenon 2 --runs 16 --seed 466
run $v L2 --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2
run $v 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551
run $v 5H --mode statistical --phenomenon 5 --runs 50 --seed 551 --pair proton,electron
run $v 6 --mode statistical --phenomenon 6 --runs 200 --seed 42
done
for tag in 1 2 L2 5Ps 5H 6; do s=SAME; for x in out.txt err.txt; do cmp -s $S/reg/old/$tag/$x $S/reg/new/$tag/$x || s=DIFF; done; echo "old vs new $tag $s"; done > $D/out/regression752.txt
echo ALLDONE

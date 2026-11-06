#!/bin/bash
# Audit 724: old = git archive 1ed9c112, new = working tree, new_off = working tree built with -DCREM_CODATA2018_CONSTANTS.  Sequential.
R=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r724
rm -rf $R/old $R/new $R/new_off $R/reg; mkdir -p $R/old; ( cd /home/teddy/Projekty/CREM-QI && git archive 1ed9c112 | tar -x -C $R/old )
cd /home/teddy/Projekty/CREM-QI && rsync -a --exclude .git --exclude audits ./ $R/new/ && rsync -a --exclude .git --exclude audits ./ $R/new_off/
( cd $R/old && make -j4 crem_interactions >/dev/null 2>&1 ); ( cd $R/new && make -j4 crem_interactions crem_interactions_validation >/dev/null 2>&1 )
( cd $R/new_off && make -j4 crem_interactions CXX="g++ -DCREM_CODATA2018_CONSTANTS" >/dev/null 2>&1 )
( cd $R/new && ./crem_interactions_validation > ../validation.txt 2>&1; echo EXIT $? >> ../validation.txt )
run(){ v=$1; tag=$2; shift 2; mkdir -p $R/reg/$v/$tag; ( cd $R/reg/$v/$tag && $R/$v/crem_interactions "$@" > out.txt 2> err.txt ); }
for v in old new_off new; do
run $v 1 --mode statistical --phenomenon 1 --runs 16 --seed 466
run $v 2 --mode statistical --phenomenon 2 --runs 16 --seed 466
run $v L2 --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2
run $v 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551
run $v 5H --mode statistical --phenomenon 5 --runs 50 --seed 551 --pair proton,electron
run $v 6 --mode statistical --phenomenon 6 --runs 200 --seed 42
done
for cmpv in new_off new; do for tag in 1 2 L2 5Ps 5H 6; do s=SAME; for x in out.txt err.txt; do cmp -s $R/reg/old/$tag/$x $R/reg/$cmpv/$tag/$x || s=DIFF; done; echo "old vs $cmpv $tag $s"; done; done > $R/regression724.txt
echo DONE

#!/bin/bash
R=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r739
rm -rf $R/old $R/new $R/reg; mkdir -p $R/old; ( cd /home/teddy/Projekty/CREM-QI && git archive 637edfd9 | tar -x -C $R/old )
cd /home/teddy/Projekty/CREM-QI && rsync -a --exclude .git --exclude audits ./ $R/new/
( cd $R/old && make -j4 crem_interactions >/dev/null 2>&1 ); ( cd $R/new && make -j4 crem_interactions crem_interactions_validation >/dev/null 2>&1 )
( cd $R/new && ./crem_interactions_validation > ../validation.txt 2>&1; echo EXIT $? >> ../validation.txt )
run(){ v=$1; tag=$2; shift 2; mkdir -p $R/reg/$v/$tag; ( cd $R/reg/$v/$tag && $R/${v%%_*}/crem_interactions "$@" > out.txt 2> err.txt ); }
for v in old new new_off; do
  if [ $v = new_off ]; then export CREM_CLASSICAL_S_TRIPLET=1; else unset CREM_CLASSICAL_S_TRIPLET; fi
run $v 1 --mode statistical --phenomenon 1 --runs 16 --seed 466
run $v 2 --mode statistical --phenomenon 2 --runs 16 --seed 466
run $v L2 --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2
run $v 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551
run $v 5H --mode statistical --phenomenon 5 --runs 50 --seed 551 --pair proton,electron
run $v 6 --mode statistical --phenomenon 6 --runs 200 --seed 42
done
unset CREM_CLASSICAL_S_TRIPLET
for cmpv in new new_off; do for tag in 1 2 L2 5Ps 5H 6; do s=SAME; for x in out.txt err.txt; do cmp -s $R/reg/old/$tag/$x $R/reg/$cmpv/$tag/$x || s=DIFF; done; echo "old vs $cmpv $tag $s"; done; done > $R/regression739.txt
echo DONE

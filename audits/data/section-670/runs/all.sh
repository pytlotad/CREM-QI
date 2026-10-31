#!/bin/bash
R=$R670
cd $R && R=$R xargs -P 4 -I{} bash -c "{}" < jobs.txt
cd /home/teddy/Projekty/CREM-QI && rsync -a --exclude .git --exclude audits ./ $R/new/
( cd $R/old && make -j4 crem_interactions >/dev/null 2>&1 ); ( cd $R/new && make -j4 crem_interactions crem_interactions_validation >/dev/null 2>&1 )
( cd $R/new && ./crem_interactions_validation > ../validation.txt 2>&1; echo EXIT $? >> ../validation.txt )
run(){ v=$1; tag=$2; shift 2; mkdir -p $R/reg/$v/$tag; ( cd $R/reg/$v/$tag && $R/$v/crem_interactions "$@" > out.txt 2> err.txt ); }
for v in old new; do
run $v 1 --mode statistical --phenomenon 1 --runs 16 --seed 466
run $v 2 --mode statistical --phenomenon 2 --runs 16 --seed 466
run $v L2 --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2
run $v 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551
run $v 5H --mode statistical --phenomenon 5 --runs 50 --seed 551 --pair proton,electron
run $v 6 --mode statistical --phenomenon 6 --runs 200 --seed 42
done
echo DONE

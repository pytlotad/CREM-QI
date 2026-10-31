#!/bin/bash
R=$R669
( cd $R/old && make -j2 crem_interactions >/dev/null 2>&1 ); ( cd $R/new && make -j2 crem_interactions >/dev/null 2>&1 )
run(){ v=$1; tag=$2; shift 2; mkdir -p $R/reg/$v/$tag; ( cd $R/reg/$v/$tag && $R/$v/crem_interactions "$@" > out.txt 2> err.txt ); }
for v in old new; do
run $v 1 --mode statistical --phenomenon 1 --runs 16 --seed 466 &
run $v 2 --mode statistical --phenomenon 2 --runs 16 --seed 466 &
run $v L2 --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2 &
run $v 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551 &
wait
run $v 5H --mode statistical --phenomenon 5 --runs 50 --seed 551 --pair proton,electron &
run $v 6 --mode statistical --phenomenon 6 --runs 200 --seed 42 &
wait
done
echo DONE

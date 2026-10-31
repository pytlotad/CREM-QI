#!/bin/bash
R=$R669
( cd $R/old && make -j2 crem_interactions >/dev/null 2>&1 ); ( cd $R/new && make -j2 crem_interactions >/dev/null 2>&1 )
run(){ v=$1; tag=$2; shift 2; mkdir -p $R/reg2/$v/$tag; ( cd $R/reg2/$v/$tag && $R/$v/crem_interactions "$@" > out.txt 2> err.txt ); }
for v in old new; do
run $v 1 --mode statistical --phenomenon 1 --runs 16 --seed 466
run $v 2 --mode statistical --phenomenon 2 --runs 16 --seed 466
run $v 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551
done
echo DONE

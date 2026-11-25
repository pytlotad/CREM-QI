#!/bin/bash
# Audit 855: CREM_DEBUG_COLLAPSE_TRANSIT (print only).  K = HEAD 1f285071 o-Ps n = 2; tree OFF o-Ps n = 2; ON p-Ps / o-Ps x n = 1 / 2 (8 runs, seed 466).
export LC_ALL=C; ulimit -v 4000000
R=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r855
rm -rf $R; mkdir -p $R/old $R/reg; ( cd /home/teddy/Projekty/CREM-QI && git archive 1f285071 | tar -x -C $R/old )
cd /home/teddy/Projekty/CREM-QI && rsync -a --exclude .git --exclude audits ./ $R/new/
( cd $R/old && make -j3 crem_interactions >/dev/null 2>&1 ) || { echo BLAD old; exit 1; }
( cd $R/new && make -j3 crem_interactions >/dev/null 2>&1 ) || { echo BLAD new; exit 1; }
echo "build $(date +%T)"
run(){ v=$1; tag=$2; shift 2; mkdir -p $R/reg/$tag; ( ulimit -v 2500000; cd $R/reg/$tag && $R/$v/crem_interactions "$@" > out.txt 2> err.txt; echo $? > rc.txt ); echo "done $tag rc $(cat $R/reg/$tag/rc.txt) $(date +%T)"; }
unset CREM_DEBUG_COLLAPSE_TRANSIT
run old K_L2 --mode statistical --phenomenon 2 --runs 8 --seed 466 --level 2
run new OFF_L2 --mode statistical --phenomenon 2 --runs 8 --seed 466 --level 2
export CREM_DEBUG_COLLAPSE_TRANSIT=1
for p in 1 2; do for n in 1 2; do run new ON_p${p}_n${n} --mode statistical --phenomenon $p --runs 8 --seed 466 --level $n; done; done
echo "ALLDONE $(date +%T)"

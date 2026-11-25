#!/bin/bash
# Audit 854: CREM_PERIAPSIS_SO_THOMAS in dipoleAwarePeriapsis.  K = HEAD ea6df0a3 set 2; tree OFF / ON sets 1, 2, L2 (32 runs, seed 466).
export LC_ALL=C; ulimit -v 4000000
R=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r854
rm -rf $R; mkdir -p $R/old $R/reg; ( cd /home/teddy/Projekty/CREM-QI && git archive ea6df0a3 | tar -x -C $R/old )
cd /home/teddy/Projekty/CREM-QI && rsync -a --exclude .git --exclude audits ./ $R/new/
( cd $R/old && make -j3 crem_interactions >/dev/null 2>&1 ) || { echo BLAD old; exit 1; }
( cd $R/new && make -j3 crem_interactions >/dev/null 2>&1 ) || { echo BLAD new; exit 1; }
echo "build $(date +%T)"
run(){ v=$1; tag=$2; shift 2; mkdir -p $R/reg/$tag; ( ulimit -v 2500000; cd $R/reg/$tag && $R/$v/crem_interactions "$@" > out.txt 2> err.txt ) && echo "done $tag $(date +%T)" || echo "FAIL $tag $(date +%T)"; }
run old K_2 --mode statistical --phenomenon 2 --runs 32 --seed 466
for m in OFF ON; do
  if [ $m = ON ]; then export CREM_PERIAPSIS_SO_THOMAS=1; else unset CREM_PERIAPSIS_SO_THOMAS; fi
  run new ${m}_1 --mode statistical --phenomenon 1 --runs 32 --seed 466
  run new ${m}_2 --mode statistical --phenomenon 2 --runs 32 --seed 466
  run new ${m}_L2 --mode statistical --phenomenon 2 --runs 32 --seed 466 --level 2
done
echo "ALLDONE $(date +%T)"

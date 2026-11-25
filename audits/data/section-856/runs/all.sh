#!/bin/bash
# Audit 856 (default change): REG old = git archive 3da85aa4, new = tree, new_legacy = tree + CREM_COLLAPSE_TRANSIT_LEGACY=1; 6 sets as 852 + exit codes;
# validation old / new; C = new --circular-start, p-Ps / o-Ps x --level 1 / 2, 8 runs, CREM_DEBUG_COLLAPSE_TRANSIT=1.
export LC_ALL=C; ulimit -v 4000000
R=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r856
rm -rf $R; mkdir -p $R/old $R/reg; ( cd /home/teddy/Projekty/CREM-QI && git archive 3da85aa4 | tar -x -C $R/old )
cd /home/teddy/Projekty/CREM-QI && rsync -a --exclude .git --exclude audits ./ $R/new/
( cd $R/old && make -j3 crem_interactions crem_interactions_validation >/dev/null 2>&1 ) || { echo BLAD old; exit 1; }
( cd $R/new && make -j3 crem_interactions crem_interactions_validation >/dev/null 2>&1 ) || { echo BLAD new; exit 1; }
echo "build $(date +%T)"
for v in old new; do ( ulimit -v 2500000; cd $R/$v && ./crem_interactions_validation > $R/validation_$v.txt 2>&1; echo EXIT $? >> $R/validation_$v.txt ); echo "validation $v $(date +%T)"; done
run(){ v=$1; tag=$2; shift 2; mkdir -p $R/reg/$v/$tag; ( ulimit -v 2500000; cd $R/reg/$v/$tag; s=$(date +%s); $R/${v%%_*}/crem_interactions "$@" > out.txt 2> err.txt; echo $? > rc.txt; echo $(( $(date +%s) - s )) > wall.txt ); }
for v in old new new_legacy; do
  if [ $v = new_legacy ]; then export CREM_COLLAPSE_TRANSIT_LEGACY=1; else unset CREM_COLLAPSE_TRANSIT_LEGACY; fi
  run $v 1 --mode statistical --phenomenon 1 --runs 16 --seed 466
  run $v 2 --mode statistical --phenomenon 2 --runs 16 --seed 466
  run $v L2 --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2
  run $v 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551
  run $v 5H --mode statistical --phenomenon 5 --runs 50 --seed 551 --pair proton,electron
  run $v 6 --mode statistical --phenomenon 6 --runs 200 --seed 42
  echo "reg $v $(date +%T)"
done
unset CREM_COLLAPSE_TRANSIT_LEGACY
export CREM_DEBUG_COLLAPSE_TRANSIT=1
for p in 1 2; do for n in 1 2; do run new C_p${p}_n${n} --mode statistical --phenomenon $p --runs 8 --seed 466 --level $n --circular-start; echo "C p$p n$n rc $(cat $R/reg/new/C_p${p}_n${n}/rc.txt) $(date +%T)"; done; done
echo "ALLDONE $(date +%T)"

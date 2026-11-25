#!/bin/bash
# Audit 860 (dynamics switch CREM_JUMP_ECCENTRIC_SHAPE, default off): REG old = git archive of the 859 commit (OLDREV), new, new_on; 6 sets as 856-859;
# validation new / new_on; D = new_on p-Ps / o-Ps --level 2 (16 runs, + ORBIT858) and --level 1 (16 runs, default budget); C = new_on --level 2 --circular-start (8 runs).
export LC_ALL=C; ulimit -v 4000000
OLDREV=${OLDREV:?}
R=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r860
rm -rf $R; mkdir -p $R/old $R/reg; ( cd /home/teddy/Projekty/CREM-QI && git archive $OLDREV | tar -x -C $R/old )
cd /home/teddy/Projekty/CREM-QI && rsync -a --exclude .git --exclude audits ./ $R/new/
( cd $R/old && make -j3 crem_interactions crem_interactions_validation >/dev/null 2>&1 ) || { echo BLAD old; exit 1; }
( cd $R/new && make -j3 crem_interactions crem_interactions_validation >/dev/null 2>&1 ) || { echo BLAD new; exit 1; }
echo "build old=$OLDREV $(date +%T)"
( ulimit -v 2500000; cd $R/new && ./crem_interactions_validation > $R/validation_new.txt 2>&1; echo EXIT $? >> $R/validation_new.txt ); echo "validation new $(date +%T)"
( ulimit -v 2500000; cd $R/new && CREM_JUMP_ECCENTRIC_SHAPE=1 ./crem_interactions_validation > $R/validation_new_on.txt 2>&1; echo EXIT $? >> $R/validation_new_on.txt ); echo "validation new_on $(date +%T)"
run(){ v=$1; tag=$2; shift 2; mkdir -p $R/reg/$v/$tag; ( ulimit -v 2500000; cd $R/reg/$v/$tag; s=$(date +%s); $R/${v%%_*}/crem_interactions "$@" > out.txt 2> err.txt; echo $? > rc.txt; echo $(( $(date +%s) - s )) > wall.txt ); }
for v in old new new_on; do
  if [ $v = new_on ]; then export CREM_JUMP_ECCENTRIC_SHAPE=1; else unset CREM_JUMP_ECCENTRIC_SHAPE; fi
  run $v 1 --mode statistical --phenomenon 1 --runs 16 --seed 466
  run $v 2 --mode statistical --phenomenon 2 --runs 16 --seed 466
  run $v L2 --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2
  run $v 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551
  run $v 5H --mode statistical --phenomenon 5 --runs 50 --seed 551 --pair proton,electron
  run $v 6 --mode statistical --phenomenon 6 --runs 200 --seed 42
  echo "reg $v $(date +%T)"
done
export CREM_JUMP_ECCENTRIC_SHAPE=1 CREM_DEBUG_COLLAPSE_TRANSIT=1
for p in 1 2; do CREM_DEBUG_TRANSIT_ORBIT=1 run new_on D_p${p}_n2 --mode statistical --phenomenon $p --runs 16 --seed 466 --level 2; echo "D p$p n2 rc $(cat $R/reg/new_on/D_p${p}_n2/rc.txt) $(date +%T)"; done
for p in 1 2; do run new_on D_p${p}_n1 --mode statistical --phenomenon $p --runs 16 --seed 466 --level 1; echo "D p$p n1 rc $(cat $R/reg/new_on/D_p${p}_n1/rc.txt) $(date +%T)"; done
for p in 1 2; do run new_on C_p${p}_n2 --mode statistical --phenomenon $p --runs 8 --seed 466 --level 2 --circular-start; echo "C p$p n2 rc $(cat $R/reg/new_on/C_p${p}_n2/rc.txt) $(date +%T)"; done
echo "ALLDONE $(date +%T)"

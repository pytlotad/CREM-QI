#!/bin/bash
# Audit 861 (860 eccentric jump shape made default; CREM_JUMP_ECCENTRIC_SHAPE_LEGACY=1 restores): REG new, new_legacy; 6 sets as 856-860;
# references = archive of audit 860 (reg/new_on, reg/new, D_p1_n2); validation new / new_legacy; D = new p-Ps --level 2 (16 runs, debug lines).
export LC_ALL=C; ulimit -v 4000000; unset CREM_JUMP_ECCENTRIC_SHAPE CREM_JUMP_ECCENTRIC_SHAPE_LEGACY CREM_DEBUG_COLLAPSE_TRANSIT CREM_DEBUG_TRANSIT_ORBIT
R=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r861
rm -rf $R; mkdir -p $R/ref860 $R/reg; tar -xJf /home/teddy/Projekty/CREM-QI/audits/data/section-860/runs/runs860.tar.xz -C $R/ref860
cd /home/teddy/Projekty/CREM-QI && rsync -a --exclude .git --exclude audits ./ $R/new/
( cd $R/new && make -j3 crem_interactions crem_interactions_validation >/dev/null 2>&1 ) || { echo BLAD new; exit 1; }
echo "build $(date +%T)"
( ulimit -v 2500000; cd $R/new && ./crem_interactions_validation > $R/validation_new.txt 2>&1; echo EXIT $? >> $R/validation_new.txt ); echo "validation new $(date +%T)"
( ulimit -v 2500000; cd $R/new && CREM_JUMP_ECCENTRIC_SHAPE_LEGACY=1 ./crem_interactions_validation > $R/validation_new_legacy.txt 2>&1; echo EXIT $? >> $R/validation_new_legacy.txt ); echo "validation new_legacy $(date +%T)"
run(){ v=$1; tag=$2; shift 2; mkdir -p $R/reg/$v/$tag; ( ulimit -v 2500000; cd $R/reg/$v/$tag; s=$(date +%s); $R/new/crem_interactions "$@" > out.txt 2> err.txt; echo $? > rc.txt; echo $(( $(date +%s) - s )) > wall.txt ); }
for v in new new_legacy; do
  if [ $v = new_legacy ]; then export CREM_JUMP_ECCENTRIC_SHAPE_LEGACY=1; else unset CREM_JUMP_ECCENTRIC_SHAPE_LEGACY; fi
  run $v 1 --mode statistical --phenomenon 1 --runs 16 --seed 466
  run $v 2 --mode statistical --phenomenon 2 --runs 16 --seed 466
  run $v L2 --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2
  run $v 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551
  run $v 5H --mode statistical --phenomenon 5 --runs 50 --seed 551 --pair proton,electron
  run $v 6 --mode statistical --phenomenon 6 --runs 200 --seed 42
  echo "reg $v $(date +%T)"
done
unset CREM_JUMP_ECCENTRIC_SHAPE_LEGACY
CREM_DEBUG_COLLAPSE_TRANSIT=1 CREM_DEBUG_TRANSIT_ORBIT=1 run new D_p1_n2 --mode statistical --phenomenon 1 --runs 16 --seed 466 --level 2
echo "D p1 n2 rc $(cat $R/reg/new/D_p1_n2/rc.txt) $(date +%T)"
echo "ALLDONE $(date +%T)"

#!/bin/bash
# Audit 858 (print-only switch CREM_DEBUG_TRANSIT_ORBIT): K = file control old (git archive becfce94) vs new without the switch,
# --phenomenon 2 --level 2 and --phenomenon 1 --level 1, 4 runs, seed 466, CREM_DEBUG_COLLAPSE_TRANSIT=1;
# P = new with the switch, p-Ps / o-Ps x --level 1 / 2; C = new with the switch, p-Ps / o-Ps --level 2 --circular-start.
export LC_ALL=C; ulimit -v 4000000
R=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r858
rm -rf $R; mkdir -p $R/old $R/reg; ( cd /home/teddy/Projekty/CREM-QI && git archive becfce94 | tar -x -C $R/old )
cd /home/teddy/Projekty/CREM-QI && rsync -a --exclude .git --exclude audits ./ $R/new/
( cd $R/old && make -j3 crem_interactions >/dev/null 2>&1 ) || { echo BLAD old; exit 1; }
( cd $R/new && make -j3 crem_interactions >/dev/null 2>&1 ) || { echo BLAD new; exit 1; }
echo "build $(date +%T)"
run(){ v=$1; tag=$2; shift 2; mkdir -p $R/reg/$v/$tag; ( ulimit -v 2500000; cd $R/reg/$v/$tag; s=$(date +%s); $R/$v/crem_interactions "$@" > out.txt 2> err.txt; echo $? > rc.txt; echo $(( $(date +%s) - s )) > wall.txt ); }
export CREM_DEBUG_COLLAPSE_TRANSIT=1
for v in old new; do
  run $v K_p2_n2 --mode statistical --phenomenon 2 --runs 4 --seed 466 --level 2
  run $v K_p1_n1 --mode statistical --phenomenon 1 --runs 4 --seed 466 --level 1
  echo "K $v $(date +%T)"
done
export CREM_DEBUG_TRANSIT_ORBIT=1
for p in 1 2; do for n in 1 2; do run new P_p${p}_n${n} --mode statistical --phenomenon $p --runs 4 --seed 466 --level $n; echo "P p$p n$n rc $(cat $R/reg/new/P_p${p}_n${n}/rc.txt) $(date +%T)"; done; done
for p in 1 2; do run new C_p${p}_n2 --mode statistical --phenomenon $p --runs 4 --seed 466 --level 2 --circular-start; echo "C p$p n2 rc $(cat $R/reg/new/C_p${p}_n2/rc.txt) $(date +%T)"; done
echo "ALLDONE $(date +%T)"

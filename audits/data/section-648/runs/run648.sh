#!/bin/bash
R=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r648; O=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r647
cd $R
{
B="CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 MODE=est QSPIN=1 BUDGET=3000"
for i in 0 1 2 3; do
 for st in "c1 1 1.0" "c2 2 1.0" "c3 3 1.0" "L1 1 0.5" "L2 2 0.75" "L3 3 0.8333333333" "n23 2.3 0.8695652174"; do set -- $st
  echo "$B PAIR=electron,positron PHEN=2 NLEV=$2 LFRAC=$3 $R/p623 $i >/dev/null 2>$R/runs/oPs_$1_$i.err"
 done
 for st in "c1 1 1.0" "c2 2 1.0" "c3 3 1.0" "n23 2.3 0.8695652174"; do set -- $st
  echo "$B PAIR=electron,positron PHEN=2 NLEV=$2 LFRAC=$3 $O/p623 $i >/dev/null 2>$R/runs/old_oPs_$1_$i.err"
 done
 for st in "c2 2 1.0" "L2 2 0.75"; do set -- $st
  echo "$B PAIR=proton,electron PHEN=1 NLEV=$2 LFRAC=$3 $R/p623 $i >/dev/null 2>$R/runs/H_$1_$i.err"
 done
 echo "$B PAIR=proton,electron PHEN=1 NLEV=2 LFRAC=1.0 $O/p623 $i >/dev/null 2>$R/runs/old_H_c2_$i.err"
done
for i in 0 1; do for st in "3d 3 0.7333333333" "3p 3 0.4"; do set -- $st
  echo "$B PAIR=electron,positron PHEN=2 NLEV=$2 LFRAC=$3 $R/p623 $i >/dev/null 2>$R/runs/oPs_$1_$i.err"; done; done
} > jobs.txt
xargs -P 6 -I{} bash -c "{}" < jobs.txt
( cd new && ./crem_interactions_validation > ../validation.txt 2>&1; echo EXIT $? >> ../validation.txt )
run(){ tag=$1; shift; mkdir -p reg/$tag; ( cd reg/$tag && $R/new/crem_interactions "$@" > out.txt 2> err.txt ); }
run 1 --mode statistical --phenomenon 1 --runs 16 --seed 466
run 2 --mode statistical --phenomenon 2 --runs 16 --seed 466
run L2 --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2
run 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551
run 5H --mode statistical --phenomenon 5 --runs 50 --seed 551 --pair proton,electron
run 6 --mode statistical --phenomenon 6 --runs 200 --seed 42
echo DONE > done.txt

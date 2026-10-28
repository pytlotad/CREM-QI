#!/bin/bash
S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r647
( cd $S/old && make crem_interactions >/dev/null 2>&1 )
cd $S
{ for st in "3d 3 0.7333333333" "3p 3 0.4" "L09 2 0.45"; do set -- $st; for i in 0 1; do echo "CREM_DEBUG_WAVE_LOCK=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 MODE=est PAIR=electron,positron PHEN=2 QSPIN=1 NLEV=$2 LFRAC=$3 BUDGET=3000 $S/p623 $i > /dev/null 2> $S/probe_$1_$i.err"; echo "CREM_WAVE_LOCK=0 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 MODE=est PAIR=electron,positron PHEN=2 QSPIN=1 NLEV=$2 LFRAC=$( [ $1 = 3d ] && echo 0.8333333333 || ([ $1 = 3p ] && echo 0.5 || echo 0.75)) BUDGET=3000 $S/p623 $i > /dev/null 2> $S/probe_${1}_ref_$i.err"; done; done; } > probes.txt
xargs -P 4 -I{} bash -c "{}" < probes.txt
( cd new && ./crem_interactions_validation > ../validation.txt 2>&1; echo EXIT $? >> ../validation.txt )
( cd new && CREM_WAVE_LOCK=0 ./crem_interactions_validation > ../validation_off.txt 2>&1; echo EXIT $? >> ../validation_off.txt )
run(){ v=$1; tag=$2; shift 2; mkdir -p reg/$v$tag; ( cd reg/$v$tag && $S/$v/crem_interactions "$@" > out.txt 2> err.txt ); }
for v in old new; do
  run $v 1 --mode statistical --phenomenon 1 --runs 16 --seed 466
  run $v 2 --mode statistical --phenomenon 2 --runs 16 --seed 466
  run $v L2 --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2
  run $v 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551
  run $v 5H --mode statistical --phenomenon 5 --runs 50 --seed 551 --pair proton,electron
  run $v 6 --mode statistical --phenomenon 6 --runs 200 --seed 42
done
echo DONE > done.txt

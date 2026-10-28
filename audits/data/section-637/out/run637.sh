#!/bin/bash
S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r637
cd $S/new && ./crem_interactions_validation > $S/validation_new.txt 2>&1; echo EXIT $? >> $S/validation_new.txt
cd $S/new && CREM_HIDDEN_J_MECHANICAL=1 ./crem_interactions_validation > $S/validation_legacy.txt 2>&1; echo EXIT $? >> $S/validation_legacy.txt
for v in old new; do for ph in 1 2; do mkdir -p $S/reg/$v$ph; ( cd $S/reg/$v$ph && $S/$v/crem_interactions --mode statistical --phenomenon $ph --runs 16 --seed 466 > out.txt 2> err.txt ); done
  mkdir -p $S/reg/${v}L2; ( cd $S/reg/${v}L2 && $S/$v/crem_interactions --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2 > out.txt 2> err.txt ); done
for j in 0 1 2 3; do ( CREM_EMISSION_LEDGER=1 MODE=est PAIR=electron,positron PHEN=2 QSPIN=1 NLEV=3 LFRAC=0.8333333333 BUDGET=3000 $S/p623 $((j*4)) $((j*4+1)) $((j*4+2)) $((j*4+3)) > $S/C_$j.out 2> $S/C_$j.err ) & done; wait
echo DONE > $S/done.txt

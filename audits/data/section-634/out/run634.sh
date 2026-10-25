#!/bin/bash
S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r634
cd /home/teddy/Projekty/CREM-QI
( cd $S/new && ./crem_interactions_validation > $S/validation.txt 2>&1; echo EXIT $? >> $S/validation.txt ) &
for ph in 1 2; do for v in old new; do mkdir -p $S/reg/$v$ph; ( cd $S/reg/$v$ph && $S/$v/crem_interactions --mode statistical --phenomenon $ph --runs 16 --seed 466 > out.txt 2> err.txt ) & done; done
wait
for v in old new; do mkdir -p $S/reg/${v}L2; ( cd $S/reg/${v}L2 && $S/$v/crem_interactions --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2 > out.txt 2> err.txt ) & done
for j in 0 1; do ( CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 MODE=est PAIR=electron,positron PHEN=2 QSPIN=1 NLEV=3 LFRAC=0.8333333333 BUDGET=3000 $S/p623 $(seq $((j*8)) $((j*8+7))) > $S/C_$j.out 2> $S/C_$j.err ) & done
wait
echo DONE > $S/done.txt

#!/bin/bash
S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r631
cd /home/teddy/Projekty/CREM-QI
( cd $S/new && ./crem_interactions_validation > $S/validation.txt 2>&1; echo EXIT $? >> $S/validation.txt ) &
for ph in 1 2; do for v in old new; do mkdir -p $S/reg/$v$ph; ( cd $S/reg/$v$ph && $S/$v/crem_interactions --mode statistical --phenomenon $ph --runs 16 --seed 466 --level 2 > out.txt 2> err.txt ) & done; done
wait
{ for ph in 1 2; do echo "== ph $ph"; cmp -s $S/reg/old$ph/out.txt $S/reg/new$ph/out.txt && echo identical || diff $S/reg/old$ph/out.txt $S/reg/new$ph/out.txt | head; done; grep -E "Enforced|verdict|EXIT" $S/validation.txt; } > $S/regression.txt
# T2 of audit 623 (id cases) with the new probe
mkdir -p $S/T2id
t2(){ cs=$1; B=$2; P0=electron,positron; PHEN=1; QSPIN=1; case $cs in ops) PHEN=2;; free) QSPIN=0;; h) P0=proton,electron;; esac
  CREM_EMISSION_LEDGER=1 MODE=est PAIR=$P0 PHEN=$PHEN QSPIN=$QSPIN NLEV=2 LFRAC=0.75 BEXT=$B BSIGN=1 BUDGET=3000 $S/p623new 0 1 2 3 > $S/T2id/${cs}_B${B}_id.out 2>$S/T2id/${cs}_B${B}_id.err; }
for cs in pps ops; do for B in 0 0.1; do t2 $cs $B & done; done; wait
for cs in free h; do for B in 0 0.1; do t2 $cs $B & done; done; wait
for f in $S/T2id/*.err; do b=$(basename $f); cmp -s $f audits/data/section-623/T2/$b && echo "T2 $b identical" || echo "T2 $b DIFFERS"; done >> $S/regression.txt
# H631c (80 seeds, 4 jobs) and LEGACY check (8 seeds, new with CREM_PHOTON_J_BASINS=1 vs old)
mkdir -p $S/C
for j in 0 1 2 3; do ( CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 MODE=est PAIR=electron,positron PHEN=2 QSPIN=1 NLEV=3 LFRAC=0.8333333333 BUDGET=3000 $S/p623new $(seq $((j*20)) $((j*20+19))) > $S/C/new_$j.out 2> $S/C/new_$j.err ) & done; wait
( CREM_PHOTON_J_BASINS=1 CREM_EMISSION_LEDGER=1 MODE=est PAIR=electron,positron PHEN=2 QSPIN=1 NLEV=3 LFRAC=0.8333333333 BUDGET=3000 $S/p623new 0 1 2 3 > $S/C/leg_a.out 2> $S/C/leg_a.err ) &
( CREM_EMISSION_LEDGER=1 MODE=est PAIR=electron,positron PHEN=2 QSPIN=1 NLEV=3 LFRAC=0.8333333333 BUDGET=3000 $S/p623old 0 1 2 3 > $S/C/old_a.out 2> $S/C/old_a.err ) &
( CREM_PHOTON_J_BASINS=1 CREM_EMISSION_LEDGER=1 MODE=est PAIR=electron,positron PHEN=2 QSPIN=1 NLEV=3 LFRAC=0.8333333333 BUDGET=3000 $S/p623new 4 5 6 7 > $S/C/leg_b.out 2> $S/C/leg_b.err ) &
( CREM_EMISSION_LEDGER=1 MODE=est PAIR=electron,positron PHEN=2 QSPIN=1 NLEV=3 LFRAC=0.8333333333 BUDGET=3000 $S/p623old 4 5 6 7 > $S/C/old_b.out 2> $S/C/old_b.err ) &
wait
for x in a b; do cmp -s $S/C/leg_$x.err $S/C/old_$x.err && echo "LEGACY $x identical" || echo "LEGACY $x DIFFERS"; done >> $S/regression.txt
echo DONE >> $S/regression.txt

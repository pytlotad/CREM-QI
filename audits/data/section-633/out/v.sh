#!/bin/bash
S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r633
cd /home/teddy/Projekty/CREM-QI
v=$1; case $v in V1) export CREM_NO_HIDDEN_ANGULAR_MOMENTUM=1;; V2) export CREM_NO_SPIN_ENERGY_EXCHANGE=1;; V3) export CREM_NO_HIDDEN_ANGULAR_MOMENTUM=1 CREM_NO_SPIN_ENERGY_EXCHANGE=1;; esac
CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 MODE=est PAIR=electron,positron PHEN=2 QSPIN=1 NLEV=3 LFRAC=0.8333333333 BUDGET=3000 $S/p623 0 3 5 > $S/$v.out 2> $S/$v.err

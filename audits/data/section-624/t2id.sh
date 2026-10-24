#!/bin/bash
S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r624
cd /home/teddy/Projekty/CREM-QI
cs=$1; B=$2
P0=electron,positron; PHEN=1; QSPIN=1
case $cs in ops) PHEN=2;; free) QSPIN=0;; h) P0=proton,electron;; esac
mkdir -p $S/T2id
CREM_EMISSION_LEDGER=1 MODE=est PAIR=$P0 PHEN=$PHEN QSPIN=$QSPIN NLEV=2 LFRAC=0.75 BEXT=$B BSIGN=1 BUDGET=3000 \
  $S/p624probe 0 1 2 3 > $S/T2id/${cs}_B${B}_id.out 2>$S/T2id/${cs}_B${B}_id.err

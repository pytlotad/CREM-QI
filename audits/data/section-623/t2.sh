#!/bin/bash
# T2: production estimator with the emission ledger; args: case B T
S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r623
cd /home/teddy/Projekty/CREM-QI
cs=$1; B=$2; T=$3
P0=electron,positron; PC=positron,electron; PHEN=1; QSPIN=1
case $cs in ops) PHEN=2;; free) QSPIN=0;; h) P0=proton,electron; PC=antiproton,positron;; esac
case $T in id) PAIR=$P0; SYM=; BS=1;; P) PAIR=$P0; SYM=P; BS=1;; C) PAIR=$PC; SYM=M; BS=-1;; CP) PAIR=$PC; SYM=PM; BS=-1;; esac
mkdir -p $S/T2
if [ -n "$SYM" ]; then export CREM_SYMMETRY_TEST=$SYM; fi
CREM_EMISSION_LEDGER=1 MODE=est PAIR=$PAIR PHEN=$PHEN QSPIN=$QSPIN NLEV=2 LFRAC=0.75 BEXT=$B BSIGN=$BS BUDGET=${BUDGET:-900} \
  $S/p623 ${IDX:-0 1 2 3} > $S/T2/${cs}_B${B}_${T}.out 2>$S/T2/${cs}_B${B}_${T}.err

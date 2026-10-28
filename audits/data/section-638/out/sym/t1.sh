#!/bin/bash
# T1: one process per (pair, field, transform); args: sys B T
S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r638/sym
cd /home/teddy/Projekty/CREM-QI
sys=$1; B=$2; T=$3
case $sys in ps) P0=electron,positron; PC=positron,electron;; h) P0=proton,electron; PC=antiproton,positron;; esac
case $T in id) PAIR=$P0; SYM=; BS=1;; P) PAIR=$P0; SYM=P; BS=1;; C) PAIR=$PC; SYM=M; BS=-1;; CP) PAIR=$PC; SYM=PM; BS=-1;; esac
mkdir -p $S/T1
if [ -n "$SYM" ]; then export CREM_SYMMETRY_TEST=$SYM; fi
MODE=traj PAIR=$PAIR BEXT=$B BSIGN=$BS FRAMES=50 $S/p623 0 1 2 3 > $S/T1/${sys}_B${B}_${T}.txt 2>$S/T1/${sys}_B${B}_${T}.err

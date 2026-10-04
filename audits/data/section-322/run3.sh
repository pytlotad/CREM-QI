#!/bin/bash
S="$1"; cd /home/teddy/Projekty/Positronium
export CREM_SPIN_TRIM=1 CREM_EMISSION_REACH=1 CREM_INITIAL_ANGULAR_MOMENTUM=0.5
CREM_LADDER_BELOW_2=1 timeout 20000 ./positronium --mode statistical --phenomenon 1 --level 2 --radiation-reaction stochastic \
    --bohr-photon-energy --ground-state-floor --runs 24 --seed 42 --crem-wallclock-budget-s 2000 > "$S/sec322/bohr.txt" 2>&1
timeout 20000 ./positronium --mode statistical --phenomenon 1 --level 2 --radiation-reaction stochastic \
    --runs 24 --seed 42 --crem-wallclock-budget-s 900 > "$S/sec322/hw.txt" 2>&1
touch "$S/sec322/DONE"

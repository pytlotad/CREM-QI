#!/bin/bash
# Audit 786 runs: probe p786 (= ../section-699/p699.cpp) built from the working tree; clock 701 with CREM_PAIR_CLOCK_SPLIT (S) and boosts 1e-3 plane (BP) / axis (BA).
# Hydrogen, antihydrogen (PHEN 2), p-Ps (PHEN 1), o-Ps (PHEN 2); seeds 0, 1; 20 radial periods, 2048 steps per orbit; H and anti-H at BEXT 0.01 T (the 683 trace needs a field), Ps without field.  3 parallel, ulimit -v 2.5 GB, LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); B=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r786/p786; O=$D/out/runs; mkdir -p $O
export B O
C="CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 MODE=est QSPIN=1 BUDGET=3000 NLEV=2 LFRAC=0.75 CREM_PAIR_QM1S_TRACE=1e-12 CREM_PAIR_CLOCK=20 CREM_MECHANICAL_STEPS_PER_ORBIT=2048 CREM_PAIR_CLOCK_SPLIT=1"
echo "start $(date +%T)"
{ for s in 0 1; do for sys in "H PAIR=proton,electron PHEN=2 BEXT=0.01" "aH PAIR=antiproton,positron PHEN=2 BEXT=0.01" "pPs PHEN=1" "oPs PHEN=2"; do set -- $sys; n=$1; shift
    echo "S_${n}_$s $C $* ./ $s"; echo "BP_${n}_$s $C CREM_PAIR_CLOCK_BOOST=1e-3,plane $* ./ $s"; echo "BA_${n}_$s $C CREM_PAIR_CLOCK_BOOST=1e-3,axis $* ./ $s"; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=${@: -1}; envs=("${@:1:$#-2}"); env "${envs[@]}" $B $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)"'
echo "ALLDONE $(date +%T)"

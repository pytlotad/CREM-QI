#!/bin/bash
# usage: probes.sh BIN OUTDIR [ENV...]  (as 734: CREM_EMISSION_LEDGER=1, deterministic threshold of p585)
B=$1; O=$2; shift 2; mkdir -p $O
run(){ tag=$1; shift; ( cd $O && env CREM_EMISSION_LEDGER=1 "$@" $B 0 > $tag.out 2> $tag.err ); }
run H2p PAIR=proton,electron NLEV=2 LFRAC=0.75 "$@"
run H3p PAIR=proton,electron NLEV=3 LFRAC=0.5 "$@"
run H4p PAIR=proton,electron NLEV=4 LFRAC=0.375 "$@"
run H3s PAIR=proton,electron NLEV=3 LFRAC=0.16666666666666666 "$@"
run H3d PAIR=proton,electron NLEV=3 LFRAC=0.8333333333333334 "$@"
run Ps2p NLEV=2 LFRAC=0.75 "$@"

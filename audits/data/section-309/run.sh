#!/bin/bash
S="$1"; T="$S/t309"
export CREM_PROBE_SPIN_QUANT=1 CREM_NO_SPIN_MAGNITUDE=1 CREM_RETARDED_BETWEEN_PHOTONS=1 CREM_PRESCRIBED_EMISSION_PATTERN=1
timeout 3000 "$T/pol" 24 42 300 > "$T/A.out" 2>&1
CREM_LS_FREE_MAGNITUDE=1 timeout 3000 "$T/pol" 24 42 300 > "$T/B.out" 2>&1
timeout 3000 "$T/s0.20/pol" 4 42 600 > "$T/A020.out" 2>&1
timeout 3000 "$T/s0.10/pol" 4 42 900 > "$T/A010.out" 2>&1
touch "$T/DONE"

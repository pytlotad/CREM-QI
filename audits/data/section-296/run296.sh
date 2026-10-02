#!/bin/bash
# R296a: przed trzecim fotonem (n > 1/2) zero odmow ceiling; po nim
# (n = 0,477 < 1/2) kazda proba emisji odmowiona przez ceiling.
S="$1"
CREM_ACTION_TRACE=1 timeout 400 "$S/stopcause" 1 300 1 \
    2> >(python3 -u "$S/stamp.py" > "$S/trace300.err") > "$S/out300.txt" &
timeout 150 "$S/stopcause" 1 80 1 > "$S/out80.txt" 2>/dev/null &
wait

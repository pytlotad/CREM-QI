#!/bin/bash
# Audit 755 runs (sequential).  LC_ALL=C.
set -u; export LC_ALL=C
D=$(cd "$(dirname "$0")/.." && pwd); O=$D/out/runs; mkdir -p $O
python3 $D/bohm755.py 0.05 1000 20000 $O/main_h005.npz
python3 $D/bohm755.py 0.025 200 2000 $O/ctrl_h0025.npz
echo done

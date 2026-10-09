#!/usr/bin/env python3
"""Audit 395 (H395b): energy change of Scattering events from CREM_CAPTURE_DUMP lines.
Usage: scatter395.py <run output>...  Prints the median |E_coul,end - E_coul,start| and its ratio to K."""
import re, sys, statistics as st
pat = re.compile(r"CAPDUMP out=Scattering K=(\S+) .*? Ecoul=(\S+)/(\S+) .*? rad=(\S+)")
for path in sys.argv[1:]:
    rows = [tuple(map(float, m.groups())) for m in map(pat.search, open(path)) if m]
    if not rows:
        print(f"== {path}: no Scattering lines"); continue
    d = [abs(c1 - c0) for K, c0, c1, rad in rows]
    rel = [abs(c1 - c0)/K for K, c0, c1, rad in rows]
    print(f"== {path}: {len(rows)} scattering events; median |dE_coul| = {st.median(d):.4g} eV, "
          f"median |dE_coul|/K = {st.median(rel):.4g}, 90% = {sorted(rel)[int(0.9*len(rel))]:.4g}")

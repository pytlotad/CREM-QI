#!/usr/bin/env python3
"""Audit 394: energy budget of the captures from CREM_CAPTURE_DUMP lines.
Usage: analyze394.py <run output>...  Energies in eV."""
import re, sys, statistics as st

pat = re.compile(r"CAPDUMP out=(\S+) K=(\S+) b=(\S+) rmin=(\S+) Ecoul=(\S+)/(\S+) Edd=(\S+)/(\S+) "
                 r"Ecd=(\S+)/(\S+) Edar=(\S+)/(\S+) Econ=(\S+)/(\S+) rad=(\S+) L=(\S+) t=(\S+)")
for path in sys.argv[1:]:
    rows = [m.groups() for m in map(pat.search, open(path)) if m]
    ev = [dict(out=r[0], K=float(r[1]), b=float(r[2]), rmin=float(r[3]),
               c0=float(r[4]), c1=float(r[5]), dd0=float(r[6]), dd1=float(r[7]),
               cd0=float(r[8]), cd1=float(r[9]), da0=float(r[10]), da1=float(r[11]),
               cn0=float(r[12]), cn1=float(r[13]), rad=float(r[14])) for r in rows]
    cap = [e for e in ev if e["out"] in ("Para-Positronium", "Ortho-Positronium")]
    print(f"== {path}: {len(ev)} events, {len(cap)} captures")
    if not cap: continue
    drop = [e["c0"] - e["c1"] for e in cap]                     # Coulomb-energy drop
    other = [(e["dd1"]+e["cd1"]+e["da1"]+e["cn1"]) - (e["dd0"]+e["cd0"]+e["da0"]+e["cn0"]) for e in cap]
    rad = [e["rad"] for e in cap]
    full_end = [e["c1"]+e["dd1"]+e["cd1"]+e["da1"]+e["cn1"] for e in cap]
    full_start = [e["c0"]+e["dd0"]+e["cd0"]+e["da0"]+e["cn0"] for e in cap]
    radfrac = [r/d for r, d in zip(rad, drop) if d > 0]
    print(f"  Coulomb drop K - E_coul,end: median {st.median(drop):.4g} eV")
    print(f"  radiated / Coulomb drop: median {st.median(radfrac):.4g}, max {max(radfrac):.4g}")
    print(f"  change of omitted terms (dd+cd+Darwin+constraint): median {st.median(other):.4g} eV")
    for name, key in (("dipole-dipole", ("dd0","dd1")), ("charge-dipole", ("cd0","cd1")),
                      ("Darwin", ("da0","da1")), ("constraint", ("cn0","cn1"))):
        d = [e[key[1]] - e[key[0]] for e in cap]
        print(f"    {name:14s} end-start: median {st.median(d):+.4g} eV, end value median {st.median([e[key[1]] for e in cap]):+.4g} eV")
    print(f"  full model energy start: median {st.median(full_start):.4g} eV; end: median {st.median(full_end):.4g} eV")
    margin = [0.05*e["K"] for e in cap]
    bound_full = sum(1 for f, m in zip(full_end, margin) if f < -m)
    print(f"  captures still bound by the FULL energy (< -0.05 K): {bound_full} of {len(cap)}")
    print(f"  r_min of captures: median {st.median([e['rmin'] for e in cap]):.4g} m")

# Audit 698 POST HOC: D_mech(t) ~ c D_sec(t) at t = 1, 2 ps -- least-squares scale c and the residual rms per run (same pattern, other rate?).
# usage: python3 post698.py <698 runs dir with mech/>
import sys, math, os
import numpy as np
from ana698 import parse, at, PAIRS, FIELDS
R = sys.argv[1]; out = ['per run: c = <D_mech, D_sec>/<D_sec, D_sec> on t = 1, 2 ps (12 states), rms residual |D_mech - c D_sec| / rms |D_mech|, c at 1 ps and 2 ps separately']
for tag, name in PAIRS:
    for B in FIELDS:
        for s in ('0', '1'):
            d = parse(f'{R}/mech/{tag}_B{B}_A2_{s}.err'); f0 = d['form']
            Dm = np.array([at(d['mech'], x) - f0 for x in (1e-12, 2e-12)]); Ds = np.array([at(d['sec'], x) - f0 for x in (1e-12, 2e-12)])
            c = float((Dm*Ds).sum()/(Ds*Ds).sum()); res = math.sqrt(((Dm - c*Ds)**2).mean())/math.sqrt((Dm**2).mean())
            c1 = float((Dm[0]*Ds[0]).sum()/(Ds[0]*Ds[0]).sum()); c2 = float((Dm[1]*Ds[1]).sum()/(Ds[1]*Ds[1]).sum())
            out.append(f'{name:12s} B {B:>6} T seed {s}: c {c:.3f} (1 ps {c1:.3f}, 2 ps {c2:.3f}) residual {res:.2f}')
text = '\n'.join(out); print(text); os.makedirs('out', exist_ok=True); open('out/post698.txt', 'w').write(text + '\n')

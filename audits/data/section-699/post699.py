# Audit 699 POST HOC: charge conjugation at zero field -- |dX| of hydrogen over antihydrogen (same seed, same field) per path, and the
# internal (B0) rates of each path over the mean field of H_2P (fs + hfs).
# usage: python3 post699.py <699 runs dir with def/>
import sys, os
import numpy as np
from ana699 import parse, hartree, FIELDS


def b0(f):
    d = parse(f); g = d['inc']['sec']; X0 = [g['L0'], g['s10'], g['s20']]
    r = {k: [np.linalg.norm(np.asarray(d['inc'][k][v]) - X0[i]) for i, v in enumerate(('L1', 's11', 's21'))] for k in ('secB0', 'mechB0')}
    hz = hartree(d['H']['fs'] + d['H']['hfs'], X0, g['dt']); r['HB0'] = [np.linalg.norm(np.asarray(hz[i]) - X0[i]) for i in range(3)]
    return r


R = sys.argv[1]; out = ['zero field (B0) |dX| ratios; vectors L^, heavy spin (p / pbar), lepton spin (e- / e+)']
for B in FIELDS:
    for s in ('0', '1'):
        h = b0(f'{R}/def/proton_B{B}_A2_{s}.err'); a = b0(f'{R}/def/antiproton_B{B}_A2_{s}.err')
        line = f'B {B:>6} T seed {s} | hydrogen / antihydrogen: ' + ' '.join(f'{k} ' + '/'.join(f'{x/y:.3f}' for x, y in zip(h[k], a[k])) for k in ('secB0', 'mechB0', 'HB0'))
        line += ' | path / Hartree: ' + ' '.join(f'{n} {k} ' + '/'.join(f'{x/y:.3f}' for x, y in zip(r[k], r['HB0'])) for n, r in (('H', h), ('Hbar', a)) for k in ('secB0', 'mechB0'))
        out.append(line)
text = '\n'.join(out); print(text); os.makedirs('out', exist_ok=True); open('out/post699.txt', 'w').write(text + '\n')
# POST HOC (2): the two transport rules of the light spin at the start state -- Born inverse cube R = L^2/(L^2 - hbar^2/4) and the
# fine-structure midpoint (528) (l + 1/2) hbar/|J_e|, J_e = L + hbar s^/2 -- against secB0 / mechB0 for L^ and the lepton spin.
import glob
out2 = ['secB0 / mechB0 against R (l + 1/2) hbar/|J_e| at the start (L^ and lepton spin), zero field']
for f in sorted(glob.glob(f'{R}/def/*_B0.01_*.err')):
    d = parse(f); g = d['inc']['sec']; Ln = g['Ln0']
    Je = np.linalg.norm(Ln*g['L0'] + 0.5*g['s20']); Rb = Ln*Ln/(Ln*Ln - 0.25); pred = Rb*1.5/Je
    obs = [np.linalg.norm(np.asarray(d['inc']['secB0'][v]) - g[v0])/np.linalg.norm(np.asarray(d['inc']['mechB0'][v]) - g[v0]) for v, v0 in (('L1', 'L0'), ('s21', 's20'))]
    out2.append(f'{os.path.basename(f)[:-4]:26s} L.s {g["L0"] @ g["s20"]:+.3f} |J_e| {Je:.4f} R {Rb:.4f} predicted {pred:.4f} observed L^ {obs[0]:.4f} lepton {obs[1]:.4f}')
text2 = '\n'.join(out2); print(text2); open('out/post699.txt', 'a').write(text2 + '\n')

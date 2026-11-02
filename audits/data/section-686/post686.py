# Audit 686, POST HOC (after ana686.py): (1) exact Hotelling p from F(3, n-3) instead of the pre-registered chi^2(3) approximation;
# (2) antihydrogen split into seeds 0-31 (already seen in 684/685) and 32-127 (new); (3) pooled per PHEN.
# usage: python3 post686.py <686 runs dir> <684 runs dir> <685 runs dir> (the ana686 arguments)
import sys, glob, math
import numpy as np
from scipy.stats import chi2, f as fdist
# ana686 runs on import (rewriting out/ana686.txt identically), so the same three directories are passed
R = sys.argv[1]
from ana686 import parse, PAIRS


def hot(rows):
    D = np.array([[r['sc'][i] - r['qm'][i] for i in range(4)] for r in rows if 'sc' in r]); n = len(D)
    mu = D.mean(0); S = np.cov(D[:, :3], rowvar=False); t2 = float(n*mu[:3] @ np.linalg.solve(S, mu[:3]))
    return n, mu, t2, float(chi2.sf(t2, 3)), float(fdist.sf(t2*(n - 3)/(3*(n - 1)), 3, n - 3))


def seed(f): return int(f.rsplit('_', 1)[1].split('.')[0])


out = ['Hotelling per cell: T^2, p chi^2(3) (pre-registered), p F exact; pooled per PHEN chi^2(6)']
for tag, name in PAIRS:
    for path in ('def', 'born'):
        perP = {'1': 0.0, '2': 0.0}
        for B in ('0.01', '0.1'):
            for P in ('1', '2'):
                files = sorted(glob.glob(f'{R}/{path}/{tag}_B{B}_A{P}_*.err'))
                n, mu, t2, pc, pf = hot([parse(x) for x in files]); perP[P] += t2
                line = f'{name:12s} {path:4s} B {B:>4} T PHEN {P} T2 {t2:6.2f} p chi2 {pc:.3g} p F {pf:.3g}'
                if tag == 'antiproton':
                    for lab, sel in (('seeds 0-31', lambda s: s < 32), ('seeds 32-127', lambda s: s >= 32)):
                        n2, mu2, t22, pc2, _ = hot([parse(x) for x in files if sel(seed(x))])
                        line += f' | {lab} n {n2} diff ' + ' '.join(f'{x:+.3f}' for x in mu2) + f' T2 {t22:5.2f} p {pc2:.3g}'
                out.append(line)
        out.append(f'{name:12s} {path:4s} pooled PHEN 1 chi^2(6) {perP["1"]:.1f} p {chi2.sf(perP["1"], 6):.3g};'
                   f' PHEN 2 chi^2(6) {perP["2"]:.1f} p {chi2.sf(perP["2"], 6):.3g}')
text = '\n'.join(out); print(text)
open('out/post686.txt', 'w').write(text + '\n')
# (4) component (A) of 684: Born rule of the photon-time state (v3) against the formation state (QM), the limit of the born path.
out2 = ['(A) = mean v3 - P_QM (Born of the photon-time state vs formation state), Hotelling chi^2(3); born path sc - v3']
for tag, name in PAIRS:
    for B in ('0.01', '0.1'):
        for P in ('1', '2'):
            rows = []
            for x in sorted(glob.glob(f'{R}/born/{tag}_B{B}_A{P}_*.err')):
                v3 = qm = sc = None
                for line in open(x):
                    t = line.split()
                    if line.startswith('PAIR QM1S684'): i = t.index('v3'); v3 = [float(y) for y in t[i+1:i+5]]
                    elif line.startswith('PAIR QM1S path'): qm = [float(y) for y in t[t.index('qm')+1:t.index('qm')+5]]; sc = [float(y) for y in t[t.index('sc')+1:t.index('sc')+5]]
                if v3 and qm: rows.append((v3, qm, sc))
            for lab, k in (('v3 - QM', lambda r: [a - b for a, b in zip(r[0], r[1])]), ('sc - v3', lambda r: [a - b for a, b in zip(r[2], r[0])])):
                D = np.array([k(r) for r in rows]); n = len(D); mu = D.mean(0)
                t2 = float(n*mu[:3] @ np.linalg.solve(np.cov(D[:, :3], rowvar=False), mu[:3]))
                out2.append(f'{name:12s} B {B:>4} T PHEN {P} {lab} n {n} ' + ' '.join(f'{x:+.4f}' for x in mu) + f'  T2 {t2:6.2f} p {chi2.sf(t2, 3):.3g}')
text2 = '\n'.join(out2); print(text2)
open('out/post686.txt', 'a').write(text2 + '\n')
# (5) born path sc - v3 split: (C) sc - v1 (the 668/670 j rule at the drawn Delta m) and v1 - v3 (draw noise; zero in expectation by the 684 identity).
out3 = ['born path: (C) sc - v1 and draw v1 - v3, Hotelling chi^2(3)']
for tag, name in PAIRS:
    for B in ('0.01', '0.1'):
        for P in ('1', '2'):
            rows = []
            for x in sorted(glob.glob(f'{R}/born/{tag}_B{B}_A{P}_*.err')):
                v = {}
                for line in open(x):
                    t = line.split()
                    if line.startswith('PAIR QM1S684'):
                        for key in ('v1', 'v3'): i = t.index(key); v[key] = [float(y) for y in t[i+1:i+5]]
                    elif line.startswith('PAIR QM1S path'): i = t.index('sc'); v['sc'] = [float(y) for y in t[i+1:i+5]]
                if len(v) == 3: rows.append(v)
            for lab, a, b in (('sc - v1', 'sc', 'v1'), ('v1 - v3', 'v1', 'v3')):
                D = np.array([[x - y for x, y in zip(r[a], r[b])] for r in rows]); n = len(D); mu = D.mean(0)
                t2 = float(n*mu[:3] @ np.linalg.solve(np.cov(D[:, :3], rowvar=False), mu[:3]))
                out3.append(f'{name:12s} B {B:>4} T PHEN {P} {lab} n {n} ' + ' '.join(f'{x:+.4f}' for x in mu) + f'  T2 {t2:6.2f} p {chi2.sf(t2, 3):.3g}')
text3 = '\n'.join(out3); print(text3)
open('out/post686.txt', 'a').write(text3 + '\n')

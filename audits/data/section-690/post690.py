# Audit 690, POST HOC (after ana690.py): split of the born path's sc - QM at 1.0329 T into (A) v3 - QM (Born rule of the photon-time state against
# the formation state) and sc - v3 (Delta m draw; no j rule at this field, so sc = v1), Hotelling chi^2(3); plus the import with the booked dE.
# usage: python3 post690.py <690 runs dir with born/ def/>
import sys, glob, math
import numpy as np
from scipy.stats import chi2
R = sys.argv[1]
UEV = 1.602176634e-25
PAIRS = (('antiproton', 'antihydrogen'), ('proton', 'hydrogen'))


def vec(t, key, n=4):
    i = t.index(key); return [float(x) for x in t[i+1:i+1+n]]


def parse(path):
    d = {}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR QM1S684'): d['v3'] = vec(t, 'v3'); d['v1'] = vec(t, 'v1')
        elif line.startswith('PAIR QM1S path'): d['sc'], d['qm'] = vec(t, 'sc'), vec(t, 'qm')
        elif line.startswith('PAIR BORNDM'): d['dE'] = float(t[t.index('dE')+1])/UEV
        elif line.startswith('PAIR eigen ledger'):
            e = {t[i]: float(t[i+1]) for i in range(3, len(t)-1, 2)}; d['imp'] = (e['B'] + e['S2P'] - e['S1S'])/UEV
    return d


def hot(D):
    D = np.array(D); n = len(D); mu = D.mean(0)
    return mu, float(n*mu[:3] @ np.linalg.solve(np.cov(D[:, :3], rowvar=False), mu[:3]))


out = ['born path at 1.0329 T: (A) v3 - QM, draw sc - v3, sc - v1 (0 without the j rule); import rms def / born / born + dE [ueV]']
for tag, name in PAIRS:
    for P in ('1', '2'):
        rows = [parse(f) for f in sorted(glob.glob(f'{R}/born/{tag}_B1.0329_A{P}_*.err'))]
        dfl = [parse(f) for f in sorted(glob.glob(f'{R}/def/{tag}_B1.0329_A{P}_*.err'))]
        for lab, a, b in (('v3 - QM', 'v3', 'qm'), ('sc - v3', 'sc', 'v3')):
            mu, t2 = hot([[x - y for x, y in zip(r[a], r[b])] for r in rows])
            out.append(f'{name:12s} PHEN {P} {lab} ' + ' '.join(f'{x:+.4f}' for x in mu) + f'  T2 {t2:6.2f} p {chi2.sf(t2, 3):.3g}')
        mx = max(max(abs(x - y) for x, y in zip(r['sc'], r['v1'])) for r in rows)
        rm = lambda xs: math.sqrt(sum(x*x for x in xs)/len(xs))
        idf = rm([r['imp'] for r in dfl]); ib = rm([r['imp'] for r in rows]); ie = rm([r['imp'] + r['dE'] for r in rows])
        out.append(f'{name:12s} PHEN {P} max |sc - v1| {mx:.1e}; import def {idf:.2f} born {ib:.2f} ({ib/idf:.2f}x) born + dE {ie:.2f} ({ie/idf:.2f}x)')
text = '\n'.join(out); print(text)
open('out/post690.txt', 'w').write(text + '\n')

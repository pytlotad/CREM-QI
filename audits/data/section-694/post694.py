# Audit 694 POST HOC: leverage of the latest photon (tau 4.45 ns) on the (A) slope, per PHEN 2 cell.
# usage: python3 post694.py <694 runs dir with poi/>
import sys, glob
import numpy as np
R = sys.argv[1]


def vec(t, key, n=4):
    i = t.index(key); return [float(x) for x in t[i+1:i+1+n]]


def parse(path):
    d = {}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR QM1S684'): d['v3'] = vec(t, 'v3')
        elif line.startswith('PAIR QM1S path'): d['qm'] = vec(t, 'qm')
        elif line.startswith('END'): d['tau'] = float.fromhex(t[t.index('tau')+1])*1e9
    return d


def fit(x, y):
    X = np.column_stack([np.ones_like(x), x]); c, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ c
    se = np.sqrt(r @ r/(len(y) - 2)*np.linalg.inv(X.T @ X)[1, 1]); return c[1], se


out = []
for c in ('antiproton_B0.01', 'antiproton_B0.1', 'proton_B0.01', 'proton_B0.1'):
    rows = [parse(f) for f in sorted(glob.glob(f'{R}/poi/{c}_A2_*.err'))]
    A = np.array([[a - b for a, b in zip(r['v3'], r['qm'])] for r in rows]); tau = np.array([r['tau'] for r in rows])
    s = int(np.argmax(np.abs(A.mean(0)))); k = tau < tau.max()
    b, se = fit(tau, A[:, s]); b2, se2 = fit(tau[k], A[k, s])
    out.append(f'{c:16s} state {s}: slope all {b:+.4f} +- {se:.4f} (t {b/se:+.2f}); without tau {tau.max():.2f} ns {b2:+.4f} +- {se2:.4f} (t {b2/se2:+.2f}); (A) at that photon {A[~k, s][0]:+.4f}')
text = '\n'.join(out); print(text); open('out/post694.txt', 'w').write(text + '\n')

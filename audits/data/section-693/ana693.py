# Audit 693 (offline): does the orbit-state term (A) = v3 - QM grow with the photon time tau (classical precession) or is it set at formation?
# usage: python3 ana693.py <686 runs dir with born/> <690 runs dir with born/>
import sys, glob, os
import numpy as np
from scipy.stats import spearmanr
R6, R9 = sys.argv[1:3]
PAIRS = (('antiproton', 'antihydrogen'), ('proton', 'hydrogen'))


def vec(t, key, n=4):
    i = t.index(key); return [float(x) for x in t[i+1:i+1+n]]


def parse(path):
    d = {}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR QM1S684'): d['v3'] = vec(t, 'v3')
        elif line.startswith('PAIR QM1S path'): d['qm'] = vec(t, 'qm')
        elif line.startswith('END'): d['tau'] = float.fromhex(t[t.index('tau')+1])
    return d


def ols(x, y):
    X = np.column_stack([np.ones_like(x), x]); c, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ c; s2 = r @ r/(len(y) - 2); se = np.sqrt(np.diag(s2*np.linalg.inv(X.T @ X)))
    return c, se


out = []; okA = okB = okC = okD = True
for B in ('0.01', '0.1', '1.0329'):
    for tag, name in PAIRS:
        for P in ('1', '2'):
            rows = [parse(f) for f in sorted(glob.glob(f'{R6 if B != "1.0329" else R9}/born/{tag}_B{B}_A{P}_*.err'))]
            rows = [r for r in rows if 'v3' in r and 'qm' in r and 'tau' in r]
            A = np.array([[a - b for a, b in zip(r['v3'], r['qm'])] for r in rows]); tau = np.array([r['tau'] for r in rows])*1e9
            s = int(np.argmax(np.abs(A.mean(0)))); (a, b), (sa, sb) = ols(tau, A[:, s])
            rho, p = spearmanr(np.linalg.norm(A, axis=1), tau); lo = tau <= np.median(tau)
            q1, q2, q3 = np.percentile(tau, (25, 50, 75))
            out.append(f'{name:12s} B {B:>6} T PHEN {P} n {len(rows)} state {s} mean {A[:, s].mean():+.4f} | a {a:+.4f} +- {sa:.4f} (t {a/sa:+.2f})'
                       f' b {b:+.4f} +- {sb:.4f} /ns (t {b/sb:+.2f}) | early {A[lo, s].mean():+.4f} late {A[~lo, s].mean():+.4f}'
                       f' | rho(|A|, tau) {rho:+.3f} p {p:.2g} | tau ns {q1:.2f} {q2:.2f} {q3:.2f} | rms |A| {np.sqrt((A**2).sum(1).mean()):.4f}')
            if B == '1.0329' and P == '2': okA &= abs(b/sb) > 3; okB &= abs(a/sa) < 2
            if B == '1.0329' and P == '1': okC &= abs(b/sb) < 3
            if P == '2': okD &= rho > 0 and p < 0.01
out.append(f'H693a (PHEN 2, 1.0329 T: |t(b)| > 3 both pairs): {okA}')
out.append(f'H693b (PHEN 2, 1.0329 T: |t(a)| < 2 both pairs): {okB}')
out.append(f'H693c (PHEN 1, 1.0329 T: |t(b)| < 3 both pairs): {okC}')
out.append(f'H693d (PHEN 2: rho(|A|, tau) > 0, p < 0.01 in all 6 cells): {okD}')
text = '\n'.join(out); print(text)
os.makedirs('out', exist_ok=True); open('out/ana693.txt', 'w').write(text + '\n')

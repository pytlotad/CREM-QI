# Audit 694: Poisson photon times (the probe default from now on) at 0.01 and 0.1 T -- does the fast dephasing of (A) = v3 - QM depend on the field,
# and the path gap sc - QM against 686 (fixed threshold).
# usage: python3 ana694.py <694 runs dir with poi/ ctl/> <686 runs dir with born/> <693 runs dir with poi/> <689 runs dir with born/>
# (control reference = 689 born, same binary as p676 with the QM1S689 line; 686 files lack that line -- corrected after the first run)
import sys, glob, math, os
import numpy as np
from scipy.stats import chi2
R, R6, R3, R89 = sys.argv[1:5]
PAIRS = (('antiproton', 'antihydrogen'), ('proton', 'hydrogen'))


def vec(t, key, n=4):
    i = t.index(key); return [float(x) for x in t[i+1:i+1+n]]


def parse(path):
    d = {}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR QM1S684'): d['v3'] = vec(t, 'v3')
        elif line.startswith('PAIR QM1S path'): d['sc'], d['qm'] = vec(t, 'sc'), vec(t, 'qm')
        elif line.startswith('END'): d['tau'] = float.fromhex(t[t.index('tau')+1])*1e9
    return d


def ols(x, y):
    X = np.column_stack([np.ones_like(x), x]); c, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ c; s2 = r @ r/(len(y) - 2); se = np.sqrt(np.diag(s2*np.linalg.inv(X.T @ X)))
    return c, se


def same(a, b): return open(a, 'rb').read() == open(b, 'rb').read()


out = []
c1 = sorted(glob.glob(f'{R}/ctl/det_*.err')); c2 = sorted(glob.glob(f'{R}/ctl/poi_*.err'))
n1 = sum(same(f, f'{R89}/born/{os.path.basename(f)[4:]}') for f in c1); n2 = sum(same(f, f'{R3}/poi/{os.path.basename(f)[4:]}') for f in c2)
out.append(f'H694a: EMIT=deterministic = 689 born {n1}/{len(c1)}; default = 693 poi {n2}/{len(c2)}: {n1 + n2 == len(c1) + len(c2) == 32}')
okB = okC = True; worst = None
for tag, name in PAIRS:
    T = 0.0
    for B in ('0.01', '0.1'):
        for P in ('2',):   # reduced scope (prereg694 uzupelnienie 1): PHEN 2 only, seeds 0-7
            rows = [r for r in (parse(f) for f in sorted(glob.glob(f'{R}/poi/{tag}_B{B}_A{P}_*.err'))) if 'v3' in r and 'sc' in r and 'tau' in r]
            det = [r for r in (parse(f) for f in sorted(glob.glob(f'{R6}/born/{tag}_B{B}_A{P}_*.err'))) if 'v3' in r and 'sc' in r]
            n = len(rows); tau = np.array([r['tau'] for r in rows])
            A = np.array([[a - b for a, b in zip(r['v3'], r['qm'])] for r in rows]); Ad = np.array([[a - b for a, b in zip(r['v3'], r['qm'])] for r in det])
            s = int(np.argmax(np.abs(A.mean(0)))); (a, b), (sa, sb) = ols(tau, A[:, s])
            e = tau < 0.15; rall = np.sqrt((A**2).sum(1).mean()); rearly = np.sqrt((A[e]**2).sum(1).mean()) if e.any() else float('nan')
            qs = np.percentile(tau, (25, 50, 75)); qi = np.digitize(tau, qs)
            G = np.array([[x - y for x, y in zip(r['sc'], r['qm'])] for r in rows]); mu = G.mean(0)
            t2 = float(n*mu[:3] @ np.linalg.solve(np.cov(G[:, :3], rowvar=False), mu[:3])); T += t2
            i = int(np.argmax(np.abs(mu))); se2 = 2*G[:, i].std(ddof=1)/math.sqrt(n)
            if P == '2': okB &= abs(b/sb) < 3   # early-photon criterion dropped (uzupelnienie 1)
            if tag == 'antiproton' and B == '0.01' and P == '2': worst = abs(mu[i])
            out.append(f'{name:12s} B {B:>4} T PHEN {P} n {n} tau mean {tau.mean():.3f} ns | (A) state {s} mean {A[:, s].mean():+.4f} (fixed 686 {Ad[:, s].mean():+.4f})'
                       f' a {a:+.4f} +- {sa:.4f} b {b:+.4f} +- {sb:.4f} /ns (t {b/sb:+.2f}) | by quartile ' + ' '.join(f'{A[qi == k, s].mean():+.4f}' for k in range(4))
                       + f' | rms |A| early (n {e.sum()}) {rearly:.4f} all {rall:.4f} (fixed {np.sqrt((Ad**2).sum(1).mean()):.4f})')
            out.append(f'   gap sc - QM ' + ' '.join(f'{x:+.4f}' for x in mu) + f' max {abs(mu[i]):.4f} (2 sigma {se2:.4f}) T2 {t2:.2f} p {chi2.sf(t2, 3):.3g}')
    p = chi2.sf(T, 6); okC &= p > 0.05; out.append(f'{name:12s} pooled gap chi^2(6) {T:.1f} p {p:.3g} (report only)')
out.append(f'H694b reduced (PHEN 2: |t(b)| < 3, both fields and pairs, n = 8): {okB}')
out.append(f'H694c not testable at n = 8 (report: gap pooled p > 0.05 both pairs {okC})')
out.append(f'H694d not testable at n = 8 (report: antihydrogen 0.01 T PHEN 2 max |gap| {worst:.4f})')
text = '\n'.join(out); print(text)
os.makedirs('out', exist_ok=True); open('out/ana694.txt', 'w').write(text + '\n')

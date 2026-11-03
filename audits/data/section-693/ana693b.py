# Audit 693 part B: Poisson photon times (EMIT=poisson, the model default) at 1.0329 T -- (A) = v3 - QM against tau, and the path gap sc - QM.
# usage: python3 ana693b.py <693 runs dir with poi/ ctl/> <690 runs dir with born/>
import sys, glob, math, os
import numpy as np
from scipy.stats import spearmanr, kstest, chi2
PAIRS = (('antiproton', 'antihydrogen'), ('proton', 'hydrogen'))


def vec(t, key, n=4):
    i = t.index(key); return [float(x) for x in t[i+1:i+1+n]]


def ols(x, y):
    X = np.column_stack([np.ones_like(x), x]); c, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ c; s2 = r @ r/(len(y) - 2); se = np.sqrt(np.diag(s2*np.linalg.inv(X.T @ X)))
    return c, se
R, R9 = sys.argv[1:3]


def parse(path):
    d = {}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR QM1S684'): d['v3'] = vec(t, 'v3')
        elif line.startswith('PAIR QM1S path'): d['sc'], d['qm'] = vec(t, 'sc'), vec(t, 'qm')
        elif line.startswith('END'): d['tau'] = float.fromhex(t[t.index('tau')+1])
    return d


out = []
ctl = sorted(glob.glob(f'{R}/ctl/*.err')); same = sum(open(f, 'rb').read() == open(f'{R9}/born/{os.path.basename(f)}', 'rb').read() for f in ctl)
out.append(f'H693e: without EMIT = 690 born byte for byte {same}/{len(ctl)}: {same == len(ctl) == 16}')
okF = okG = okH = okI = True; pooled = {}
for tag, name in PAIRS:
    T = 0.0
    for P in ('1', '2'):
        rows = [parse(f) for f in sorted(glob.glob(f'{R}/poi/{tag}_B1.0329_A{P}_*.err'))]
        rows = [r for r in rows if 'v3' in r and 'qm' in r and 'tau' in r and 'sc' in r]; n = len(rows)
        A = np.array([[a - b for a, b in zip(r['v3'], r['qm'])] for r in rows]); tau = np.array([r['tau'] for r in rows])*1e9
        s = int(np.argmax(np.abs(A.mean(0)))); (a, b), (sa, sb) = ols(tau, A[:, s])
        rho, p = spearmanr(np.linalg.norm(A, axis=1), tau)
        ks = kstest(tau, 'expon', args=(0, tau.mean())).pvalue; zt = (tau.mean() - 1.6748)/(tau.std(ddof=1)/math.sqrt(n))
        okF &= ks > 0.01 and abs(zt) < 3
        qs = np.percentile(tau, (25, 50, 75)); qi = np.digitize(tau, qs)
        quart = ' '.join(f'{A[qi == k, s].mean():+.4f}' for k in range(4))
        G = np.array([[x - y for x, y in zip(r['sc'], r['qm'])] for r in rows]); mu = G.mean(0)
        t2 = float(n*mu[:3] @ np.linalg.solve(np.cov(G[:, :3], rowvar=False), mu[:3])); T += t2
        i = int(np.argmax(np.abs(mu))); se2 = 2*G[:, i].std(ddof=1)/math.sqrt(n)
        out.append(f'{name:12s} PHEN {P} n {n} tau mean {tau.mean():.4f} ns (z vs 1.6748 {zt:+.2f}) KS p {ks:.2g} quartile edges ' + ' '.join(f'{x:.3f}' for x in qs))
        out.append(f'   (A) state {s} mean {A[:, s].mean():+.4f} | a {a:+.4f} +- {sa:.4f} (t {a/sa:+.2f}) b {b:+.4f} +- {sb:.4f} /ns (t {b/sb:+.2f})'
                   f' | rho(|A|, tau) {rho:+.3f} p {p:.2g} | (A)_s by tau quartile {quart} | rms |A| {np.sqrt((A**2).sum(1).mean()):.4f}')
        out.append(f'   gap sc - QM ' + ' '.join(f'{x:+.4f}' for x in mu) + f' max {abs(mu[i]):.4f} (2 sigma {se2:.4f}) T2 {t2:.2f} p {chi2.sf(t2, 3):.3g}')
        if P == '2': okG &= abs(b/sb) > 3; okH &= abs(a/sa) < 2
        else: okI &= abs(b/sb) < 3
    pooled[name] = (T, chi2.sf(T, 6)); out.append(f'{name:12s} pooled gap chi^2(6) {T:.1f} p {chi2.sf(T, 6):.3g}')
out.append(f'H693f (tau exponential: KS p > 0.01 and mean within 3 SE of 1.6748 ns, every cell): {okF}')
out.append(f'H693g (PHEN 2: |t(b)| > 3 both pairs): {okG}')
out.append(f'H693h (PHEN 2: |t(a)| < 2 both pairs): {okH}')
out.append(f'H693i (PHEN 1: |t(b)| < 3 both pairs): {okI}')
out.append(f'H693j (gap pooled p > 0.05 both pairs): {all(v[1] > 0.05 for v in pooled.values())}')
text = '\n'.join(out); print(text)
os.makedirs('out', exist_ok=True); open('out/ana693b.txt', 'w').write(text + '\n')

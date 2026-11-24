# Audit 844: engine spin part (Z / ZD / ZDC) vs BO'C c^-2, orbit scaling s.  Written before the runs.
import os, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(D, 'out', 'runs')
HEV, C = 4.135667696e-15, 299792458.0
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:+.4f} MHz = {f*HEV*1e6:+.5f} ueV' + (f' (|lambda| {C/abs(f)*100:.4g} cm)' if f else '')
def gauss(N):
    x, w = np.polynomial.legendre.leggauss(N); u = 0.75 + 0.25*x[::-1]; return 1/u, 0.25*w[::-1]
NODES = dict(zip([f'u5{c}' for c in 'abcde'], zip(*gauss(5))))
nu_of = lambda Om: sum(w*L**2*Om[k] for k, (L, w) in NODES.items())/(2*math.pi)
SS = [1, 2, 4, 8, 16, 32]; FORMS = ('Z', 'ZD', 'ZDC', 'BOC')
SYS = {'ops': 'o-Ps (PHEN=2)', 'pps': 'p-Ps (PHEN=1)', 'H0': 'wodor F = 0 (PHEN=1)', 'H1': 'wodor F = 1 (PHEN=2)'}
def rd(f):
    w = open(f).read().split(); return {q: float(w[w.index(q) + 1]) for q in FORMS}
P841 = {'ops': (-0.0745e6, 0.0850e6), 'H0': (-0.5259e6, 0.0753e6), 'H1': (0.5259e6, 0.0762e6)}
okA = okK = True
for sk, lab in SYS.items():
    nu = {s: {q: nu_of({nd: rd(f'{R}/{sk}_{nd}_s{s}.txt')[q] for nd in NODES}) for q in FORMS} for s in SS}
    nuN2 = {s: {q: nu_of({nd: rd(f'{R}/{sk}_{nd}_s{s}_N2.txt')[q] for nd in NODES}) for q in FORMS} for s in (1, 32)}
    pr(f'{lab}:')
    pr(f'  s = 1: nu_BOC {units(nu[1]["BOC"])}')
    for q in ('Z', 'ZD', 'ZDC'):
        pr(f'    {q:3s}: nu_silnik {units(nu[1][q])};  reszta alpha^6 {units(nu[1][q] - nu[1]["BOC"])}')
    for q in ('Z', 'ZD', 'ZDC'):
        d = np.array([nu[s][q]/nu[s]['BOC'] - 1 for s in SS]); sv = np.array(SS, float)
        A = np.vstack([np.ones_like(sv), 1/sv, sv**-1.5, sv**-2]).T; c0, c1, c15, c2 = np.linalg.lstsq(A, d, rcond=None)[0]
        ok = abs(c0) <= 1e-6 and abs(d[-1]) <= abs(d[0])/10; okA &= ok
        pr(f'    {q:3s}: delta(s) ' + ' '.join(f'{x:+.3e}' for x in d) + f'  s*delta ' + ' '.join(f'{x*s:+.3e}' for x, s in zip(d, SS))
           + f';  c0 {c0:+.2e}, c1 {c1:+.3e}  -> {"TAK" if ok else "NIE"}')
    n2 = max(abs(nuN2[s][q]/nu[s][q] - 1) for s in (1, 32) for q in FORMS if nu[s][q]); okK &= n2 <= 1e-6
    msg = f'  N x 2: max wzgl. {n2:.1e}'
    if sk in P841:
        a, b = nu[1]['ZD'] - nu[1]['Z'], nu[1]['ZDC'] - nu[1]['ZD']; okK &= abs(a - P841[sk][0]) <= 200 and abs(b - P841[sk][1]) <= 200
        msg += f';  ZD - Z {a/1e6:+.4f} (841 {P841[sk][0]/1e6:+.4f}), ZDC - ZD {b/1e6:+.4f} (843 {P841[sk][1]/1e6:+.4f}) MHz'
    pr(msg)
pr(f"H844a (GLOWNA: |c0| <= 1e-6 i |delta(32)| <= |delta(1)|/10 we wszystkich ukladach i postaciach): {'TAK' if okA else 'NIE'}")
pr(f"H844k (kontrola 841 / 843 do 0,0002 MHz; N x 2 <= 1e-6): {'TAK' if okK else 'NIE'}")

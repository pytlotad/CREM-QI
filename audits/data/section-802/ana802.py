# Audit 802: secular averages on Sommerfeld-orbit nodes (S) vs Kepler (K) and the per-orbit mechanics of 799/800.  usage: python3 ana802.py
# 799 / 800 / 801 runs unpacked into $R799, $R800, $R801 (../section-79x|80x/runs/runs*.tar.xz).
import os, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
R799 = os.environ.get('R799', SP + '/r799x/runs'); R800 = os.environ.get('R800', SP + '/r800x/runs'); R801 = os.environ.get('R801', SP + '/r801x/runs')
HEV, C, ALPHA = 4.135667696e-15, 299792458.0, 7.2973525693e-3; FS = 10969.045e6
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.2f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
TORI = [(3, 1), (3, 1.5), (3, 2), (3, 2.5), (3, 3), (2, 1), (2, 1.5), (2, 2)]
ECC = [(3, 1), (2, 1), (3, 1.5), (2, 1.5), (3, 2)]
TYPES = ('proton', 'antiproton'); SEEDS = {1: (4, 5, 6, 7, 8, 9, 10, 11, 12), 0: (4, 5, 6)}
lab = lambda L: '%g' % L
nm = lambda t, k, s: f'{t}_n{k[0]}L{lab(k[1])}_{s}.err'
def ref(F, t, k, s):
    if F == 1: return f'{R799}/A/{nm(t, k, s)}' if s in (4, 5, 6) else f'{R800}/B/{nm(t, k, s)}'
    return f'{R799}/A0/{nm(t, k, s)}' if s == 5 else f'{R800}/B0/{nm(t, k, s)}'
new = lambda F, t, k, s: f'{O}/S/F{F}_{nm(t, k, s)}'
def om(p):   # Omega_obr as in ana799 / ana800 / ana801
    L = []; S = []; per = None
    for ln in open(p):
        if not ln.startswith('PAIR TORQUE799 mechB0'): continue
        w = ln.split()
        if w[3] == 'apo':
            i = w.index('L'); L.append(np.array([float(x) for x in w[i + 1:i + 4]])); i = w.index('s2'); S.append(np.array([float(x) for x in w[i + 1:i + 4]]))
        elif w[3] == 'sum': per = float(w[w.index('period') + 1])
    d = [np.linalg.norm(S[j] - S[j - 1])/np.linalg.norm(np.cross(L[j - 1], S[j - 1])) for j in range(1, len(S))]
    return float(np.mean(d))/(2*math.pi*per)
def v3(w, tag): i = w.index(tag); return np.array([float(x) for x in w[i + 1:i + 4]])
def secline(p):
    for ln in open(p):
        if ln.startswith('PAIR QM1S699 secB0'): return ln.split()
def sec801(p):
    out = {}
    for ln in open(p):
        if ln.startswith('PAIR SEC801'): w = ln.split(); out[w[2]] = (int(w[4]), v3(w, 'second'))
    return out
def rate(omega, Lh, s): s = s/np.linalg.norm(s); Lh = Lh/np.linalg.norm(Lh); return np.linalg.norm(np.cross(omega, s))/np.linalg.norm(np.cross(Lh, s))/(2*math.pi)
est = lambda k: 1.5*ALPHA**2/k[1]**2 - ALPHA**2/(2*k[0]**2)
K = os.path.join(O, 'K', 'proton_n3L1_5.err')
okK = [ln for ln in open(K) if not ln.startswith('PAIR SEC801 S')] == list(open(f'{R801}/K/proton_n3L1_5.err'))
pr(f"H802k (K bez linii SEC801 S = 801 K/proton_n3L1_5.err bajt w bajt): {'TAK' if okK else 'NIE'}")
W = {}; okV = True
for F, seeds in SEEDS.items():
    for s in seeds:
        for t in TYPES:
            for k in TORI:
                wn = secline(new(F, t, k, s)); e = sec801(new(F, t, k, s)); okV &= all(e[x][0] == 1 for x in 'KRS')
                L0, s20 = v3(wn, 'L0'), v3(wn, 's20')
                W[(F, t, k, s)] = dict(K=rate(e['K'][1], L0, s20), S=rate(e['S'][1], L0, s20), obr=om(ref(F, t, k, s)))
pr(f"linie SEC801 K, R, S valid we wszystkich biegach: {'TAK' if okV else 'NIE'}")
pm = lambda F, k, s, f: 0.5*(f(W[(F, 'proton', k, s)]) + f(W[(F, 'antiproton', k, s)]))
okE = okS = okM = okB = True; B2 = []
for F, seeds in SEEDS.items():
    pr(f'=== F = {F} (PHEN={2 if F else 1}) ===')
    pr('  torus  | S/K - 1 (sr.)  szacunek  iloraz | p(S)/p(obr) - 1 (zakres ziaren, sr.) | p(K)/p(obr) - 1 (sr.)')
    for k in TORI:
        sk = np.mean([W[(F, t, k, s)]['S']/W[(F, t, k, s)]['K'] - 1 for s in seeds for t in TYPES]); q = sk/est(k)
        okE &= abs(q - 1) <= 0.3
        xs = [pm(F, k, s, lambda x: x['S'])/pm(F, k, s, lambda x: x['obr']) - 1 for s in seeds]
        xk = [pm(F, k, s, lambda x: x['K'])/pm(F, k, s, lambda x: x['obr']) - 1 for s in seeds]
        okS &= all(abs(x) <= 2e-5 for x in xs)
        if k in ECC: okM &= abs(np.mean(xs)) < 0.3*abs(np.mean(xk))
        pr(f'  n{k[0]}L{lab(k[1]):<4} | {sk:+.3e} {est(k):+.3e} {q:.3f} | {min(xs):+.2e} .. {max(xs):+.2e} (sr. {np.mean(xs):+.2e}) | {np.mean(xk):+.2e}')
    for n, l, st in ((2, 1, '2p'), (3, 1, '3p'), (3, 2, '3d')):
        q = []
        for s in seeds:
            od = pm(F, (n, l), s, lambda x: l*x['obr']) - pm(F, (n, l + 1), s, lambda x: (l + 1)*x['obr'])
            born = pm(F, (n, l + 0.5), s, lambda x: x['S'])*(l + 0.5)**2/(l*(l + 1)); q.append(born/od - 1)
            if st == '2p': B2.append(born)
        okB &= all(abs(x) <= 3e-4 for x in q)
        pr(f'  H802b {st}: Born x p(Omega_S(l + 1/2)) / p(Omega_rozn) - 1 = {min(q):+.2e} .. {max(q):+.2e} (sr. {np.mean(q):+.2e})')
pr(f"H802e (S/K - 1 = szacunek do 30 %): {'TAK' if okE else 'NIE'}")
pr(f"H802s (p(S)/p(obr) w 2e-5): {'TAK' if okS else 'NIE'}")
pr(f"H802m (mimosrodowe: |reszta S| < 0,3 |reszta K|): {'TAK' if okM else 'NIE'}")
pr(f"H802b (Born x S = regula mechaniki w 3e-4): {'TAK' if okB else 'NIE'}")
m = np.mean(B2); pr(f'RAPORT 2P z wezlow Sommerfelda (9/8 x p(Omega_S) na n2L1.5, {len(B2)} par): {units(m)} = {m/FS:.5f} pomiaru ({(m - FS)/1e6:+.2f} MHz); zakres {min(B2)/1e6:.2f}-{max(B2)/1e6:.2f} MHz')
pr(f'  o-Ps przez czynnik z wodoru (2P modelu / pomiar x 0,9998..0,9999): {min(B2)/FS*0.9998:.5f}-{max(B2)/FS*0.9999:.5f} multipletu (Hagena i in., PRL 71, 2887 (1993)); p-Ps bez czesci spinowej')
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)).')

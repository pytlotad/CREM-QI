# Audit 806: convergence order of the residual R = E_mech - E_S of the difference rule in the step (512 / 1024 / 2048 per orbit, seed 5).
# usage: python3 ana806.py.  512: 799 A / A0 ($R799) + S from 802 ($R802); 1024: 805 V2 ($R805); 2048: out/runs/V3.
import os, sys, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); sys.path.insert(0, D)
from per_orbit import om
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
R799 = os.environ.get('R799', SP + '/r799x/runs'); R802 = os.environ.get('R802', SP + '/r802x/runs'); R805 = os.environ.get('R805', SP + '/r805x/runs')
HEV, C = 4.135667696e-15, 299792458.0; FS = 10969.045e6
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.2f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm = {f/FS:.6f} pomiaru'
TORI = [(3, 1), (3, 1.5), (3, 2), (3, 2.5), (3, 3), (2, 1), (2, 1.5), (2, 2)]
RULES = ((2, 1, '2p'), (3, 1, '3p'), (3, 2, '3d')); TYPES = ('proton', 'antiproton'); lab = lambda L: '%g' % L
H = (512, 1024, 2048)
nm = lambda t, k: f'{t}_n{k[0]}L{lab(k[1])}_5.err'
def paths(h, F, t, k):
    if h == 512: return (f'{R799}/{"A" if F == 1 else "A0"}/{nm(t, k)}', f'{R802}/S/F{F}_{nm(t, k)}')
    p = f'{R805}/V2/F{F}_{nm(t, k)}' if h == 1024 else f'{O}/V3/F{F}_{nm(t, k)}'; return (p, p)
def v3(w, tag): i = w.index(tag); return np.array([float(x) for x in w[i + 1:i + 4]])
def sec(p):
    e = {}; wn = None
    for ln in open(p):
        if ln.startswith('PAIR SEC801'): w = ln.split(); e[w[2]] = (int(w[4]), v3(w, 'second'))
        elif ln.startswith('PAIR QM1S699 secB0') and wn is None: wn = ln.split()
    L0, s20 = v3(wn, 'L0'), v3(wn, 's20'); s = s20/np.linalg.norm(s20); Lh = L0/np.linalg.norm(L0)
    r = lambda o: np.linalg.norm(np.cross(o, s))/np.linalg.norm(np.cross(Lh, s))/(2*math.pi)
    return r(e['S'][1]), all(e[x][0] == 1 for x in 'KS')
W = {}; okV = True; rej = 0; napo = []
for h in H:
    for F in (1, 0):
        for t in TYPES:
            for k in TORI:
                pm_, ps_ = paths(h, F, t, k)
                if h == 2048 and not any(ln.startswith('PAIR TORQUE799 mechB0 sum') for ln in open(pm_)): okV = False; continue
                w, r, _, nT = om(pm_); s, v = sec(ps_)
                if h == 2048: rej += r; okV &= v; napo.append(nT)
                W[(h, F, t, k)] = dict(obr=w, S=s)
okV &= rej == 0 and len(napo) == 32 and min(napo) >= 195
pr(f"H806v (32 pliki 2048: linia mechB0 sum, SEC801 K / S valid, 0 odrzuconych, >= 195 apocentrow): {'TAK' if okV else 'NIE'} (odrzuconych {rej}, apocentrow min {min(napo) if napo else 0})")
pm = lambda h, F, k, x: 0.5*(W[(h, F, 'proton', k)][x] + W[(h, F, 'antiproton', k)][x])
def E(h, F, n, l, x):
    od = l*pm(h, F, (n, l), x) - (l + 1)*pm(h, F, (n, l + 1), x)
    return od/(pm(h, F, (n, l + 0.5), x)*(l + 0.5)**2/(l*(l + 1))) - 1, od
def extrap(y):   # y at 512, 1024, 2048 -> (monotone, p, limit)
    d1, d2 = y[0] - y[1], y[1] - y[2]; mono = d1*d2 > 0 and abs(d2) < abs(d1)
    if not mono: return False, float('nan'), float('nan')
    p = math.log2(d1/d2); return True, p, y[2] - d2/(2**p - 1)
res = {}
for F in (1, 0):
    pr(f'=== F = {F} (PHEN={F + 1}), ziarno 5, pary (wodor + antywodor)/2 ===')
    pr('  regula | R(512)     R(1024)    R(2048)    | E_mech(512) E_mech(1024) E_mech(2048) | E_S        | monot. | p     | R0')
    okM = okP = okR = True
    for n, l, st in RULES:
        R = [E(h, F, n, l, 'obr')[0] - E(h, F, n, l, 'S')[0] for h in H]; Em = [E(h, F, n, l, 'obr')[0] for h in H]; eS = E(2048, F, n, l, 'S')[0]
        mono, p, R0 = extrap(R); okM &= mono; okP &= mono and 1.5 <= p <= 4.5; okR &= abs(R[2]) <= 5e-6
        res[(F, st)] = (mono, p, R0)
        pr(f'  {st}     | {R[0]:+.3e} {R[1]:+.3e} {R[2]:+.3e} | {Em[0]:+.3e}  {Em[1]:+.3e}   {Em[2]:+.3e}   | {eS:+.3e} | {"tak" if mono else "nie"}    | {p:5.2f} | {R0:+.2e}')
    if okM and okP: okR &= all(abs(res[(F, st)][2]) <= 3e-6 for _, _, st in RULES)
    pr(f"  H806m F = {F} (d1, d2 tego samego znaku, |d2| < |d1|): {'TAK' if okM else 'NIE'}")
    pr(f"  H806p F = {F} (p w [1,5; 4,5]): {'TAK' if okP else 'NIE'}")
    pr(f"  H806r F = {F} (|R(2048)| <= 5e-6; |R0| <= 3e-6 gdy H806m i H806p): {'TAK' if okR else 'NIE'}")
    t22 = [pm(2048, F, (2, 2), 'obr')/pm(2048, F, (2, 2), 'S') - 1]
    pr(f"  H806t F = {F} (|n2L2 p(mech)/p(S) - 1| przy 2048 <= 3e-6): {'TAK' if abs(t22[0]) <= 3e-6 else 'NIE'} ({t22[0]:+.2e})")
    pr('  reszty torusow p(mech)/p(S) - 1 (pary; 512 / 1024 / 2048): ' + '; '.join(f'n{k[0]}L{lab(k[1])} ' + ' / '.join('%+.2e' % (pm(h, F, k, 'obr')/pm(h, F, k, 'S') - 1) for h in H) for k in TORI))
    P = [E(h, F, 2, 1, 'obr')[1] for h in H]; eS = E(2048, F, 2, 1, 'S')[0]
    for h, x in zip(H, P): pr(f'  2P {h:4d} krokow: regula roznicowa {units(x)}; bez czlonu E_S {units(x/(1 + eS))}')
    mono, p, P0 = extrap(P)
    if mono and 1.5 <= p <= 4.5: pr(f'  2P granica (p = {p:.2f}): regula {units(P0)}; bez E_S {units(P0/(1 + eS))}'); lim = (P0, P0)
    else:
        a = [P[2] - (P[1] - P[2])/(2**q - 1) for q in (2, 4)]; lim = (min(a), max(a))
        pr(f'  2P granica (p = 2 .. 4, rzad obserwowany nieuzyteczny): regula {units(lim[0])} .. {units(lim[1])}; bez E_S {units(lim[0]/(1 + eS))} .. {units(lim[1]/(1 + eS))}')
    if F == 1: L1 = (lim, eS)
(lo, hi), eS = L1
pr(f'  o-Ps przez czynnik z wodoru (granica, F = 1): regula {lo/FS*0.9998:.5f}-{hi/FS*0.9999:.5f}, bez E_S {lo/(1 + eS)/FS*0.9998:.5f}-{hi/(1 + eS)/FS*0.9999:.5f} multipletu (Hagena i in., PRL 71, 2887 (1993)); p-Ps bez czesci spinowej')
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)).')

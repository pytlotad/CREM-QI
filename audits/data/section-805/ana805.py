# Audit 805: ablations of the 1-3e-5 residual R = E_mech - E_S of the difference rule (seed 5).  usage: python3 ana805.py
# V0: 799 A / A0 ($R799) + S from 802 ($R802); V1 (CREM_NO_DIPOLE_FORCE=1), V2 (1024 steps): out/runs/V1, V2.
import os, sys, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); sys.path.insert(0, D)
from per_orbit import om
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
R799 = os.environ.get('R799', SP + '/r799x/runs'); R802 = os.environ.get('R802', SP + '/r802x/runs')
HEV, C = 4.135667696e-15, 299792458.0; FS = 10969.045e6
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.2f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
TORI = [(3, 1), (3, 1.5), (3, 2), (3, 2.5), (3, 3), (2, 1), (2, 1.5), (2, 2)]
RULES = ((2, 1, '2p'), (3, 1, '3p'), (3, 2, '3d')); TYPES = ('proton', 'antiproton'); lab = lambda L: '%g' % L
nm = lambda t, k: f'{t}_n{k[0]}L{lab(k[1])}_5.err'
def paths(V, F, t, k):
    if V == 'V0': return (f'{R799}/{"A" if F == 1 else "A0"}/{nm(t, k)}', f'{R802}/S/F{F}_{nm(t, k)}')
    p = f'{O}/{V}/F{F}_{nm(t, k)}'; return (p, p)
def v3(w, tag): i = w.index(tag); return np.array([float(x) for x in w[i + 1:i + 4]])
def sec(p):
    e = {}; wn = None
    for ln in open(p):
        if ln.startswith('PAIR SEC801'): w = ln.split(); e[w[2]] = (int(w[4]), v3(w, 'second'))
        elif ln.startswith('PAIR QM1S699 secB0') and wn is None: wn = ln.split()
    L0, s20 = v3(wn, 'L0'), v3(wn, 's20'); s = s20/np.linalg.norm(s20); Lh = L0/np.linalg.norm(L0)
    r = lambda o: np.linalg.norm(np.cross(o, s))/np.linalg.norm(np.cross(Lh, s))/(2*math.pi)
    return r(e['S'][1]), all(e[x][0] == 1 for x in 'KS')
W = {}; okV = True; rej = 0
for V in ('V0', 'V1', 'V2'):
    for F in (1, 0):
        for t in TYPES:
            for k in TORI:
                pm_, ps_ = paths(V, F, t, k)
                if V != 'V0' and not any(ln.startswith('PAIR TORQUE799 mechB0 sum') for ln in open(pm_)): okV = False; continue
                w, r, _, _ = om(pm_); s, v = sec(ps_); rej += r if V != 'V0' else 0; okV &= v if V != 'V0' else True
                W[(V, F, t, k)] = dict(obr=w, S=s)
okV &= rej == 0
pr(f"H805v (nowe pliki: linia mechB0 sum, SEC801 K / S valid, 0 odrzuconych odcinkow): {'TAK' if okV else 'NIE'} (odrzuconych {rej})")
pm = lambda V, F, k, x: 0.5*(W[(V, F, 'proton', k)][x] + W[(V, F, 'antiproton', k)][x])
def E(V, F, n, l, x):
    od = l*pm(V, F, (n, l), x) - (l + 1)*pm(V, F, (n, l + 1), x)
    return od/(pm(V, F, (n, l + 0.5), x)*(l + 0.5)**2/(l*(l + 1))) - 1, od
verdict = {}
for F in (1, 0):
    pr(f'=== F = {F} (PHEN={F + 1}), ziarno 5, pary (wodor + antywodor)/2 ===')
    pr('  regula | wariant | E_mech     E_S        R = E_mech - E_S | R_V/R_V0')
    R = {}
    for n, l, st in RULES:
        for V in ('V0', 'V1', 'V2'):
            em, _ = E(V, F, n, l, 'obr'); es, _ = E(V, F, n, l, 'S'); R[(V, st)] = em - es
            pr(f'  {st}     | {V}      | {em:+.3e} {es:+.3e} {em - es:+.3e}        | {(em - es)/R[("V0", st)]:+.3f}')
    d_expl = all(abs(R[('V1', st)]) <= 0.3*abs(R[('V0', st)]) for _, _, st in RULES)
    d_none = all(abs(R[('V1', st)] - R[('V0', st)]) <= 0.3*abs(R[('V0', st)]) for _, _, st in RULES)
    z_none = all(abs(R[('V2', st)] - R[('V0', st)]) <= 0.3*abs(R[('V0', st)]) for _, _, st in RULES)
    verdict[F] = (d_expl, d_none, z_none)
    pr(f"  H805d F = {F}: {'tlumaczy' if d_expl else ('nie tlumaczy' if d_none else 'czesciowo')}")
    pr(f"  H805z F = {F}: {'TAK (krok nie jest zrodlem)' if z_none else 'NIE'}")
    for st in ('2p', '3p', '3d'):
        d = R[('V2', st)] - R[('V0', st)]
        pr(f'    Richardson {st}: R(h -> 0) przy rzedzie 2 {R[("V2", st)] + d/3:+.3e}, przy rzedzie 4 {R[("V2", st)] + d/15:+.3e}')
    pr('  reszty torusow p(mech)/p(S) - 1 (pary): ' + '; '.join(f'n{k[0]}L{lab(k[1])} ' + ' / '.join('%+.2e' % (pm(V, F, k, 'obr')/pm(V, F, k, 'S') - 1) for V in ('V0', 'V1', 'V2')) for k in TORI) + '  (V0 / V1 / V2)')
    pr('  S wariantu / S V0 - 1: ' + '; '.join(f'n{k[0]}L{lab(k[1])} ' + ' / '.join('%+.1e' % (pm(V, F, k, 'S')/pm('V0', F, k, 'S') - 1) for V in ('V1', 'V2')) for k in TORI) + '  (V1 / V2)')
    for V in ('V0', 'V1', 'V2'):
        em, od = E(V, F, 2, 1, 'obr'); es, _ = E(V, F, 2, 1, 'S'); c = od/(1 + es)
        pr(f'  2P {V}: regula roznicowa {units(od)} = {od/FS:.6f} pomiaru; bez czlonu E_S {units(c)} = {c/FS:.6f} pomiaru')
        if F == 1 and V == 'V0': P = (od, c)
pr(f'  o-Ps przez czynnik z wodoru (V0, F = 1): {P[0]/FS*0.9998:.5f}-{P[0]/FS*0.9999:.5f} (regula), {P[1]/FS*0.9998:.5f}-{P[1]/FS*0.9999:.5f} (bez E_S) multipletu (Hagena i in., PRL 71, 2887 (1993)); p-Ps bez czesci spinowej')
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)).')

# Audit 807: difference rule D vs Born rule B for hydrogen 2P3/2 - 2P1/2 at the step limit, against measurement + explicit QED import.
# usage: python3 ana807.py.  2048 steps: seeds 4, 6, 7, 8 (F = 1), 4, 6 (F = 0) out/runs/V3; seed 5 from 806 ($R806).
# Step correction from seed 5: 512 (799 A / A0, $R799), 1024 (805 V2, $R805), 2048 (806 V3).
import os, sys, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); sys.path.insert(0, D)
from per_orbit import om
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
R799 = os.environ.get('R799', SP + '/r799x/runs'); R805 = os.environ.get('R805', SP + '/r805x/runs'); R806 = os.environ.get('R806', SP + '/r806x/runs')
HEV, C = 4.135667696e-15, 299792458.0
NU_NAT, S_REF = 10969.045e6, math.hypot(0.015e6, 0.020e6)
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.3f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm = {f/NU_NAT:.6f} pomiaru'
TORI = [(2, 1), (2, 1.5), (2, 2)]; TYPES = ('proton', 'antiproton'); lab = lambda L: '%g' % L
SEEDS = {1: (4, 5, 6, 7, 8), 0: (4, 5, 6)}
nm = lambda t, k, s: f'{t}_n{k[0]}L{lab(k[1])}_{s}.err'
def path(h, F, t, k, s):
    if h == 512: return f'{R799}/{"A" if F == 1 else "A0"}/{nm(t, k, s)}'
    if h == 1024: return f'{R805}/V2/F{F}_{nm(t, k, s)}'
    return f'{R806 if s == 5 else O}/V3/F{F}_{nm(t, k, s)}'
okV = True; rej = 0; napo = []; W = {}
def get(h, F, s):
    w = {}
    for t in TYPES:
        for k in TORI:
            p = path(h, F, t, k, s); new = h == 2048 and s != 5
            global okV, rej
            if new and not any(ln.startswith('PAIR TORQUE799 mechB0 sum') for ln in open(p)): okV = False; continue
            o, r, _, nT = om(p)
            if new: rej += r; napo.append(nT)
            w[(t, k)] = o
    pm = lambda k: 0.5*(w[('proton', k)] + w[('antiproton', k)])
    return {'D': 1*pm((2, 1)) - 2*pm((2, 2)), 'B': 1.125*pm((2, 1.5))}
for F in (1, 0):
    for s in SEEDS[F]: W[(2048, F, s)] = get(2048, F, s)
    for h in (512, 1024): W[(h, F, 5)] = get(h, F, 5)
okV &= rej == 0 and len(napo) == 36 and min(napo) >= 195
pr(f"H807v (36 nowych plikow: linia mechB0 sum, 0 odrzuconych, >= 195 apocentrow): {'TAK' if okV else 'NIE'} (odrzuconych {rej}, apocentrow min {min(napo) if napo else 0})")
pr(f'Cel nu_nat = 10 969,045 MHz (pomiar), import QED 0,000 +- 0,020 MHz, sigma_ref = {S_REF/1e6:.3f} MHz')
def limit(y):
    d1, d2 = y[0] - y[1], y[1] - y[2]
    if d1*d2 > 0 and abs(d2) < abs(d1):
        p = math.log2(d1/d2)
        if 1.5 <= p <= 4.5: return y[2] - d2/(2**p - 1), f'p = {p:.2f}', y[2] - d2/(2**p - 1), y[2] - d2/(2**p - 1)
    a = [y[2] - d2/(2**q - 1) for q in (2, 4)]
    return 0.5*(a[0] + a[1]), 'p = 2..4 (srodek)', min(a), max(a)
verdicts = {}; res = {}
for F in (1, 0):
    pr(f'=== F = {F} (PHEN={F + 1}); pary (wodor + antywodor)/2; ziarna {SEEDS[F]} przy 2048 krokach ===')
    st = {}
    for X in ('D', 'B'):
        y5 = [W[(h, F, 5)][X] for h in (512, 1024, 2048)]; L, how, lo, hi = limit(y5); c = L - y5[2]
        v = np.array([W[(2048, F, s)][X] for s in SEEDS[F]]); m = v.mean(); ss = v.std(ddof=1)/math.sqrt(len(v))
        sh = max(abs(c), abs(hi - lo)/2, 0.02e6); sig = math.sqrt(ss**2 + sh**2 + S_REF**2); nu = m + c; z = (nu - NU_NAT)/sig
        st[X] = (nu, sig, z)
        pr(f'  regula {X}: ziarna przy 2048 ' + ', '.join(f'{s}: {x/1e6:.3f}' for s, x in zip(SEEDS[F], v)) + ' MHz')
        pr(f'    ziarno 5: 512 {y5[0]/1e6:.3f}, 1024 {y5[1]/1e6:.3f}, 2048 {y5[2]/1e6:.3f}, granica {L/1e6:.3f} MHz ({how}); c = {c/1e6:+.3f} MHz')
        pr(f'    srednia 2048 {m/1e6:.3f} MHz, sigma_s {ss/1e6:.3f}, sigma_h {sh/1e6:.3f}, sigma {sig/1e6:.3f} MHz')
        pr(f'    nu_{X} = {units(nu)}; nu - nu_nat = {(nu - NU_NAT)/1e6:+.3f} MHz = {z:+.1f} sigma')
    zD, zB = abs(st['D'][2]), abs(st['B'][2])
    pr(f"  H807D F = {F} (|nu_D - nu_nat| <= 3 sigma): {'TAK' if zD <= 3 else 'NIE'}")
    pr(f"  H807B F = {F} (|nu_B - nu_nat| <= 3 sigma): {'TAK' if zB <= 3 else 'NIE'}")
    if zD <= 3 and zB > 5: v_ = 'D wygrywa'
    elif zB <= 3 and zD > 5: v_ = 'B wygrywa'
    elif zD > 3 and zB > 3: v_ = 'obie wykluczone'
    else: v_ = 'nierozstrzygniete'
    verdicts[F] = v_; res[F] = st
    pr(f'  rozstrzygniecie F = {F}: {v_}')
    pr(f"  nu_D/nu_B - 1 = {st['D'][0]/st['B'][0] - 1:+.3e}")
fin = verdicts[1] if verdicts[1] == verdicts[0] else 'nierozstrzygniete (F = 0 i F = 1 rozne)'
pr(f'ROZSTRZYGNIECIE (to samo dla F = 0 i F = 1 wymagane): {fin}')
for X in ('D', 'B'):
    r = res[1][X][0]/NU_NAT
    pr(f'  o-Ps przez czynnik z wodoru (regula {X}, F = 1): {r*0.9998:.5f}-{r*0.9999:.5f} multipletu (Hagena i in., PRL 71, 2887 (1993)); p-Ps bez czesci spinowej')
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)).')
pr('Import QED: Blumer i in., arXiv:2412.19580, Tab. 1; Horbatsch i Hessels, PRA 93, 022513 (2016).')

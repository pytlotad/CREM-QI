# Audit 798: 699 window ending at apoapsis (CREM_PAIR_TORQUE_APO=0.5) -- does the 1/N quantum of 797 vanish?  usage: python3 ana798.py
# 797 references unpacked from ../section-797/runs/runs797.tar.xz into $R797.
import sys, os, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
sys.path.insert(0, os.path.join(D, '..', 'section-699')); import ana699
R797 = os.environ.get('R797', '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r797x/runs')
HEV, C = 4.135667696e-15, 299792458.0; FS = 10969.045e6
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.2f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
TORI = [(3, 1), (3, 1.5), (3, 2), (3, 2.5), (3, 3), (2, 1), (2, 1.5), (2, 2)]
TYPES = ('proton', 'antiproton'); SEEDS = '456'
lab = lambda L: '%g' % L
f = lambda d, t, k, s: f'{d}/{t}_n{k[0]}L{lab(k[1])}_{s}.err'
def omega(path):
    m = ana699.parse(path)['inc']['mechB0']
    return np.linalg.norm(m['s21'] - m['s20'])/(m['dt']*np.linalg.norm(np.cross(m['L0'], m['s20'])))/(2*math.pi)
def rule(W, tag, out):
    ok_d = ok_m = True
    for n, l, st in ((2, 1, '2p'), (3, 1, '3p'), (3, 2, '3d')):
        od = l*W[(n, l)] - (l + 1)*W[(n, l + 1)]; q = od/W[(n, l + 0.5)]/((l + 0.5)**2/(l*(l + 1))); ok_d &= abs(q - 1) <= 0.01
        out.setdefault(st, []).append(q)
        if st == '2p': ok_m &= abs(od/FS - 1) <= 0.01; out.setdefault('2P', []).append(od)
        pr(f'    {tag} {st}: {units(od)} | {q:.5f} x regula' + (f' | {od/FS:.5f} pomiaru' if st == '2p' else ''))
    return ok_d, ok_m
okK = open(f'{O}/K/proton_n3L1_5.err').read() == open(f'{R797}/G/proton_n3L1_5.err').read()
pr(f"H798k (bez przelacznika = 797 G proton n3L1 ziarno 5 bajt w bajt): {'TAK' if okK else 'NIE'}")
# H798q
okQ = True; cells = []
for k in ((3, 1), (2, 1)):
    for s in SEEDS:
        for t in TYPES:
            q = (omega(f(f'{O}/AK', t, k, s))/omega(f(f'{O}/A', t, k, s)) - 1)*200; okQ &= abs(q) <= 0.3; cells.append(f'n{k[0]}L1 {t[0]}{s} {q:+.3f}')
pr('  stary start / nowy start - 1, x 200 (okno w apocentrum): ' + ', '.join(cells))
pr(f"H798q (|.| <= 0,3 dla 12 par; 797 do 2,06): {'TAK' if okQ else 'NIE'}")
okZ = True
for t in TYPES:
    q = omega(f(f'{O}/AZ', t, (3, 1), 5))/omega(f(f'{O}/A', t, (3, 1), 5)) - 1; okZ &= abs(q) <= 1e-3
    pr(f'  AZ {t} n3L1 ziarno 5: Omega(2048)/Omega(512) - 1 = {q:+.3e}')
pr(f"H798z (|2048/512 - 1| <= 1e-3): {'TAK' if okZ else 'NIE'}")
# A: rule
okP = okD = okM = okG = True; out = {}; ALL = {}; W = {}
pr('=== A (F = 1, okno w apocentrum) ===')
for s in SEEDS:
    for t in TYPES:
        W[(t, s)] = {k: omega(f(f'{O}/A', t, k, s)) for k in TORI}; W7 = {k: omega(f(f'{R797}/G', t, k, s)) for k in TORI}
        for k in TORI: ALL.setdefault(k, []).append(W[(t, s)][k])
        K = {k: W[(t, s)][k]*k[0]**3*k[1]**2 for k in TORI}; sp = max(K.values())/min(K.values()) - 1; okP &= sp <= 0.01
        pr(f'  {t} ziarno {s}: rozrzut Omega n^3 L^2 {sp*100:.3f} % | A / 797 - 1, x 200: ' + ', '.join(f'n{k[0]}L{lab(k[1])} {(W[(t, s)][k]/W7[k] - 1)*200:+.2f}' for k in TORI))
        d, m = rule(W[(t, s)], f'{t} {s}', out); okD &= d; okM &= m
    g = W[('proton', s)][(3, 1)]/W[('antiproton', s)][(3, 1)] - 1; okG &= abs(g) <= 0.005
    g7 = omega(f(f'{R797}/G', 'proton', (3, 1), s))/omega(f(f'{R797}/G', 'antiproton', (3, 1), s)) - 1
    pr(f'  luka n3L1 ziarno {s}: proton / antyproton - 1 = {g*100:+.3f} % (797: {g7*100:+.3f} %)')
pr(f"H798p (rozrzut <= 1 %): {'TAK' if okP else 'NIE'}")
pr(f"H798d (2p {min(out['2p']):.5f}-{max(out['2p']):.5f}, 3p {min(out['3p']):.5f}-{max(out['3p']):.5f}, 3d {min(out['3d']):.5f}-{max(out['3d']):.5f} w 1 %): {'TAK' if okD else 'NIE'}")
pr(f"H798m (2P {min(out['2P'])/1e6:.2f}-{max(out['2P'])/1e6:.2f} MHz = {min(out['2P'])/FS:.4f}-{max(out['2P'])/FS:.4f} pomiaru, w 1 %): {'TAK' if okM else 'NIE'}")
pr(f"H798g (luka n3L1 <= 0,5 %): {'TAK' if okG else 'NIE'}")
pr('  rozrzut miedzy 6 biegami (max/min - 1) na torusie, A wobec 797 G:')
for k in TORI:
    a = ALL[k]; b = [omega(f(f'{R797}/G', t, k, s)) for s in SEEDS for t in TYPES]
    pr(f'    n{k[0]}L{lab(k[1])}: A {max(a)/min(a) - 1:.2e} (x200 {(max(a)/min(a) - 1)*200:.2f}) | 797 {max(b)/min(b) - 1:.2e} (x200 {(max(b)/min(b) - 1)*200:.2f})')
M = {k: np.mean(v) for k, v in ALL.items()}
pr('  srednie 6 biegow:'); rule(M, 'srednie', {})
# A0
okS = True; out0 = {}; okD0 = okM0 = True
pr('=== A0 (F = 0, PHEN=1, okno w apocentrum, ziarno 5) ===')
for t in TYPES:
    W0 = {k: omega(f(f'{O}/A0', t, k, 5)) for k in TORI}
    qs = {k: W0[k]/W[(t, '5')][k] - 1 for k in TORI}; okS &= all(abs(q) <= 1e-3 for q in qs.values())
    pr(f'  {t}: Omega(F = 0)/Omega(F = 1) - 1: ' + ', '.join(f'n{k[0]}L{lab(k[1])} {q:+.2e}' for k, q in qs.items()))
    d, m = rule(W0, f'{t} F = 0', out0); okD0 &= d; okM0 &= m
pr(f"H798s (|F0/F1 - 1| <= 1e-3, wszystkie torusy): {'TAK' if okS else 'NIE'}")
pr(f"H798s0 (F = 0: 2p, 3p, 3d i 2P w 1 %): {'TAK' if okD0 and okM0 else 'NIE'}")
R21 = [x*9/8 for x in out['2p']]
pr(f'  o-Ps: 8/9 x R x 0,9998..0,9999 = {8/9*min(R21)*0.9998:.5f}-{8/9*max(R21)*0.9999:.5f} multipletu (Hagena i in., PRL 71, 2887 (1993)); p-Ps bez czesci spinowej')
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)).')

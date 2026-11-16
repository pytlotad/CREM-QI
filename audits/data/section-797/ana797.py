# Audit 797: spin-orbit difference rule WITH moment forces on the new start (795 + 796); repeat of 791; F = 0 with forces.  usage: python3 ana797.py
# 791 G references unpacked from ../section-791/runs/runs791.tar.xz into $R791.
import sys, os, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
sys.path.insert(0, os.path.join(D, '..', 'section-699')); import ana699
R791 = os.environ.get('R791', '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r791x/runs')
HEV, C = 4.135667696e-15, 299792458.0; FS = 10969.045e6
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.2f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
TORI = [(3, 1), (3, 1.5), (3, 2), (3, 2.5), (3, 3), (2, 1), (2, 1.5), (2, 2)]
lab = lambda L: '%g' % L
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
okK = open(f'{O}/K/proton_n3L1_5.err').read() == open(f'{R791}/G/proton_n3L1_5.err').read()
pr(f"H797k (CREM_KEPLER_START=1 = 791 G proton n3L1 ziarno 5 bajt w bajt): {'TAK' if okK else 'NIE'}")
okP = okD = okM = okG = True; out = {}; ALL = {}; Wall = {}
pr('=== G (F = 1, sily od momentow, nowy start) ===')
for s in '456':
    for t in ('proton', 'antiproton'):
        W = {k: omega(f'{O}/G/{t}_n{k[0]}L{lab(k[1])}_{s}.err') for k in TORI}; Wall[(t, s)] = W
        W1 = {k: omega(f'{R791}/G/{t}_n{k[0]}L{lab(k[1])}_{s}.err') for k in TORI}
        for k in TORI: ALL.setdefault(k, []).append(W[k])
        K = {k: W[k]*k[0]**3*k[1]**2 for k in TORI}; sp = max(K.values())/min(K.values()) - 1; okP &= sp <= 0.01
        pr(f'  {t} ziarno {s}: rozrzut Omega n^3 L^2 {sp*100:.3f} % | wobec 791 (stary start): ' + ', '.join(f'n{k[0]}L{lab(k[1])} {W[k]/W1[k] - 1:+.2e}' for k in TORI))
        d, m = rule(W, f'{t} {s}', out); okD &= d; okM &= m
    g = Wall[('proton', s)][(3, 1)]/Wall[('antiproton', s)][(3, 1)] - 1
    g1 = omega(f'{R791}/G/proton_n3L1_{s}.err')/omega(f'{R791}/G/antiproton_n3L1_{s}.err') - 1; okG &= abs(g) <= 0.005
    pr(f'  luka n3L1 ziarno {s}: proton / antyproton - 1 = {g*100:+.3f} % (791: {g1*100:+.3f} %)')
pr(f"H797p (rozrzut <= 1 %): {'TAK' if okP else 'NIE'}")
pr(f"H797d (2p {min(out['2p']):.5f}-{max(out['2p']):.5f}, 3p {min(out['3p']):.5f}-{max(out['3p']):.5f}, 3d {min(out['3d']):.5f}-{max(out['3d']):.5f} w 1 %): {'TAK' if okD else 'NIE'}")
pr(f"H797m (2P {min(out['2P'])/1e6:.2f}-{max(out['2P'])/1e6:.2f} MHz = {min(out['2P'])/FS:.4f}-{max(out['2P'])/FS:.4f} pomiaru, w 1 %): {'TAK' if okM else 'NIE'}")
pr(f"H797g (luka n3L1 <= 0,5 %): {'TAK' if okG else 'NIE'}")
M = {k: np.mean(v) for k, v in ALL.items()}; pr('  srednie 6 biegow:'); rule(M, 'srednie', {})
okZ = True
for t in ('proton', 'antiproton'):
    q = omega(f'{O}/Z/{t}_n3L1_5.err')/Wall[(t, '5')][(3, 1)] - 1; okZ &= abs(q) <= 1e-3
    pr(f'  Z {t} n3L1 ziarno 5: Omega(2048)/Omega(512) - 1 = {q:+.2e}')
pr(f"H797z (|2048/512 - 1| <= 1e-3): {'TAK' if okZ else 'NIE'}")
pr('=== G0 (F = 0, PHEN=1, sily od momentow, ziarno 5) ===')
okS = True; out0 = {}
for t in ('proton', 'antiproton'):
    W0 = {k: omega(f'{O}/G0/{t}_n{k[0]}L{lab(k[1])}_5.err') for k in TORI}
    pr(f'  {t}: Omega(F = 0)/Omega(F = 1) - 1: ' + ', '.join(f'n{k[0]}L{lab(k[1])} {W0[k]/Wall[(t, "5")][k] - 1:+.2e}' for k in TORI))
    d, m = rule(W0, f'{t} F = 0', out0); okS &= d and m
pr(f"H797s0 (F = 0: 2p, 3p, 3d i 2P w 1 %): {'TAK' if okS else 'NIE'}")
R21 = [x*9/8 for x in out['2p']]
pr(f'  o-Ps: 8/9 x R x 0,9998..0,9999 = {8/9*min(R21)*0.9998:.5f}-{8/9*max(R21)*0.9999:.5f} multipletu (Hagena i in., PRL 71, 2887 (1993)); p-Ps bez czesci spinowej')
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)).')

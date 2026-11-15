# Audit 791: difference rule for the spin-orbit rate (prereg repetition of 790) -- converged step (512), new seeds 4-6, tori n = 2 (L = 1, 1.5, 2) and
# n = 3 (L = 1 ... 3) with CREM_PAIR_TORQUE_ANYL; F(L) = L Omega(L), Omega_diff(n, l) = F(l) - F(l+1) against (l+1/2)^2/(l(l+1)) Omega(l+1/2) and the measured 2P.
# usage: python3 ana791.py   (runs in out/runs/{F,Z,G}; 790 reference file in $R790F, default scratchpad copy)
import sys, os, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
sys.path.insert(0, os.path.join(D, '..', 'section-699')); import ana699
R790F = os.environ.get('R790F', '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r790F/runs/proton_n2L1.5_2.err')
HEV, C = 4.135667696e-15, 299792458.0; FS = 10969.045e6
pr = lambda *a: print(*a, flush=True)
TORI = [(3, 1), (3, 1.5), (3, 2), (3, 2.5), (3, 3), (2, 1), (2, 1.5), (2, 2)]
def omega(path):
    m = ana699.parse(path)['inc']['mechB0']
    return np.linalg.norm(m['s21'] - m['s20'])/(m['dt']*np.linalg.norm(np.cross(m['L0'], m['s20'])))/(2*math.pi), m['Ln0']
lab = lambda L: '%g' % L
units = lambda f: f'{f/1e6:.2f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
ref = open(R790F).read()
okF = all(open(f'{O}/F/proton_n2L1.5_2_{w}.err').read() == ref for w in ('off', 'on'))
pr(f"H791f (kontrola pliku: bez przelacznika i z nim = 790 bajt w bajt): {'TAK' if okF else 'NIE'}")
okZ = True
for L in (2, 2.5, 3):
    a = omega(f'{O}/G/proton_n3L{lab(L)}_4.err')[0]; b = omega(f'{O}/Z/proton_n3L{lab(L)}_4.err')[0]; okZ &= abs(a/b - 1) <= 1e-3
    pr(f'  n 3 L {L}: 512 {a/1e6:.3f} MHz, 2048 {b/1e6:.3f} MHz -> {a/b - 1:+.2e}')
pr(f"H791z (zbieznosc 512 / 2048 w 1e-3): {'TAK' if okZ else 'NIE'}")
okP = okD = okM = True; R21 = []
for t in ('proton', 'antiproton'):
    for s in ('4', '5', '6'):
        pr(f'=== {t} ziarno {s} ===')
        W = {}
        for n, L in TORI:
            om, ln = omega(f'{O}/G/{t}_n{n}L{lab(L)}_{s}.err'); W[(n, L)] = (om, ln)
            pr(f'  n {n} L {L:3}: Ln0 {ln:.6f} | Omega/2pi {units(om)} | Omega n^3 L^2 {om*n**3*ln**2/1e6:.3f} MHz')
        k = [om*n**3*ln**2 for (n, L), (om, ln) in W.items()]; sp = max(k)/min(k) - 1; okP &= sp <= 0.01
        pr(f'  H791p rozrzut Omega n^3 L^2: {sp*100:.3f} %')
        for n, l, st in ((2, 1, '2p'), (3, 1, '3p'), (3, 2, '3d')):
            F = lambda L: W[(n, L)][1]*W[(n, L)][0]
            od = F(l) - F(l + 1); oc = W[(n, l + 0.5)][0]; q = od/oc; qt = (l + 0.5)**2/(l*(l + 1)); okD &= abs(q/qt - 1) <= 0.01
            rf = FS*16/(n**3*l*(l + 1))
            pr(f'  {st}: Omega_rozn/2pi {units(od)} | Omega(l+1/2) {oc/1e6:.2f} MHz | stosunek {q:.5f} wobec {qt:.5f} -> {q/qt:.5f}'
               f' | wobec skalowania pomiaru 2P {rf/1e6:.2f} MHz -> {od/rf:.5f} ({(od - rf)/1e6:+.2f} MHz)')
            if st == '2p': okM &= abs(od/FS - 1) <= 0.01; R21.append(q)
pr(f"H791p (Omega n^3 L^2 stale w 1 %): {'TAK' if okP else 'NIE'}")
pr(f"H791d (regula roznicowa = (l+1/2)^2/(l(l+1)) w 1 %, 2p 3p 3d): {'TAK' if okD else 'NIE'}")
pr(f"H791m (wodor 2P: Omega_rozn = 10 969,045(15) MHz w 1 %): {'TAK' if okM else 'NIE'}")
pr('=== RAPORT: Ps 2^3P_J (o-Ps: 8/9 x R x zegar 789; p-Ps bez czesci spinowej) ===')
pr(f'  R (2p, 6 przebiegow) {min(R21):.5f}-{max(R21):.5f}; 8/9 x R x 0,9998..0,9999 = {8/9*min(R21)*0.9998:.5f}-{8/9*max(R21)*0.9999:.5f}')
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)); Ps 2^3P_J: Hagena i in., PRL 71, 2887 (1993).')

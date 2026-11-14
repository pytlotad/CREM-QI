# Audit 790: difference rule for the spin-orbit rate on the mechanical path -- Omega(n, L) from the mechB0 lines (CREM_PAIR_TORQUE, field zeroed),
# F(L) = L Omega(L), Omega_diff(n, l) = F(l) - F(l+1) against (l+1/2)^2/(l(l+1)) x Omega(l+1/2) and the measured 2P3/2 - 2P1/2 interval of hydrogen.
# usage: python3 ana790.py   (runs in out/runs; 789 B runs for the control in $R789B, default scratchpad copy)
import sys, os, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
sys.path.insert(0, os.path.join(D, '..', 'section-699')); import ana699
R789B = os.environ.get('R789B', '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r789B/runs/B')
HEV, C = 4.135667696e-15, 299792458.0; FS = 10969.045e6
pr = lambda *a: print(*a, flush=True)
TORI = {2: (1, 1.5, 2), 3: (1, 1.5, 2, 2.5, 3)}
def omega(path):   # lepton spin precession rate [Hz] and L/hbar from the mechB0 line
    m = ana699.parse(path)['inc']['mechB0']
    return np.linalg.norm(m['s21'] - m['s20'])/(m['dt']*np.linalg.norm(np.cross(m['L0'], m['s20'])))/(2*math.pi), m['Ln0']
lab = lambda L: ('%g' % L)
units = lambda f: f'{f/1e6:.2f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
pr('=== H790k: kontrola wobec 789 B (n = 2, L = 1,5) ===')
okK = True
for t in ('proton', 'antiproton'):
    for s in ('2', '3'):
        a = [l for l in open(f'{O}/{t}_n2L1.5_{s}.err') if l.startswith('PAIR QM1S699 mechB0')]
        b = [l for l in open(f'{R789B}/{t}_{s}.err') if l.startswith('PAIR QM1S699 mechB0')]
        same = len(a) == 1 and a == b; okK &= same; pr(f'  {t} ziarno {s}: mechB0 bajt w bajt {same}')
pr(f"H790k: {'TAK' if okK else 'NIE'}")
okP = okD = okM = True; R21 = []; res = {}
for t in ('proton', 'antiproton'):
    for s in ('2', '3'):
        pr(f'=== {t} ziarno {s} ===')
        W = {}
        for n, Ls in TORI.items():
            for L in Ls:
                try: om, ln = omega(f'{O}/{t}_n{n}L{lab(L)}_{s}.err')
                except KeyError: pr(f'  n {n} L {L:3}: BRAK linii mechB0 (bramka l = 1 w crem_collapse.hpp:3420) -- nieoceniane'); continue
                W[(n, L)] = (om, ln)
                pr(f'  n {n} L {L:3}: Ln0 {ln:.6f} | Omega/2pi {units(om)} | Omega n^3 L^2 {om*n**3*ln**2/1e6:.3f} MHz')
        k = [om*n**3*ln**2 for (n, L), (om, ln) in W.items()]; sp = max(k)/min(k) - 1; okP &= sp <= 0.01
        pr(f'  (H790p na {len(W)} z 8 torusow)')
        pr(f'  H790p rozrzut Omega n^3 L^2: {sp*100:.3f} %')
        for n, l, st in ((2, 1, '2p'), (3, 1, '3p'), (3, 2, '3d')):
            if not all((n, L) in W for L in (l, l + 0.5, l + 1)): pr(f'  {st}: NIEOCENIANE (brak torusa)'); continue
            F = lambda L: W[(n, L)][1]*W[(n, L)][0]
            od = F(l) - F(l + 1); oc = W[(n, l + 0.5)][0]; q = od/oc; qt = (l + 0.5)**2/(l*(l + 1)); okD &= abs(q/qt - 1) <= 0.01
            ref = FS*16/(n**3*l*(l + 1)); res[(t, s, st)] = od
            pr(f'  {st}: Omega_rozn/2pi {units(od)} | Omega(l+1/2) {oc/1e6:.2f} MHz | stosunek {q:.5f} wobec (l+1/2)^2/(l(l+1)) = {qt:.5f} -> {q/qt:.5f}'
               f' | wobec skalowania pomiaru 2P {ref/1e6:.2f} MHz -> {od/ref:.5f}')
            if st == '2p': okM &= abs(od/FS - 1) <= 0.01; R21.append(q)
pr(f"H790p (Omega n^3 L^2 stale w 1 %): {'TAK' if okP else 'NIE'}")
pr(f"H790d (regula roznicowa = (l+1/2)^2/(l(l+1)) w 1 %; ocenialny tylko 2p, 3p i 3d bez torusow): {'TAK' if okD else 'NIE'}")
pr(f"H790m (wodor 2P: Omega_rozn = 10 969,045(15) MHz w 1 %): {'TAK' if okM else 'NIE'}")
pr('=== RAPORT: Ps 2^3P_J (zegar 789 x czynnik R z wodoru) ===')
for R in R21:
    for s, c in (('2', 0.9998), ('3', 0.9999)):
        pr(f'  R = {R:.5f}, zegar 789 ziarno {s} ({c}): o-Ps czesc spinowa / multiplet ze stalymi zmierzonymi = 8/9 x R x {c} = {8/9*R*c:.5f}')
pr('  p-Ps: brak czesci spinowej (S = 0), regula go nie dotyczy.')
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)); Ps 2^3P_J: Hagena i in., PRL 71, 2887 (1993).')

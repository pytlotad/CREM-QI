# Audit 795: production fix of the periapsis start (osculatingPeriapsisState on the relativistic torus; CREM_KEPLER_START=1 restores pre-795).
# usage: python3 ana795.py   (regression / validation copied from the scratchpad into out/; 794 references unpacked into $R794)
import sys, os, math, re
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
sys.path.insert(0, os.path.join(D, '..', 'section-699')); import ana699
R794 = os.environ.get('R794', '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r794x/runs')
HEV, C = 4.135667696e-15, 299792458.0; FS = 10969.045e6
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.2f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
TORI = [(3, 1), (3, 1.5), (3, 2), (3, 2.5), (3, 3), (2, 1), (2, 1.5), (2, 2)]
lab = lambda L: '%g' % L
def omega(path):
    m = ana699.parse(path)['inc']['mechB0']
    return np.linalg.norm(m['s21'] - m['s20'])/(m['dt']*np.linalg.norm(np.cross(m['L0'], m['s20'])))/(2*math.pi)
def relstart(path):
    for ln in open(path):
        if ln.startswith('PAIR TORQUE794'):
            w = ln.split(); return {w[i]: float(w[i + 1]) for i in range(3, len(w) - 1, 2)}
# --- H795o, H795v
reg = open(os.path.join(D, 'out', 'regression795.txt')).read().split('\n')
okO = all(l.endswith('SAME') for l in reg if l.startswith('old vs new_off'))
pr('\n'.join(l for l in reg if l)); pr(f"H795o (old = new_off bajt w bajt, 6 zestawow): {'TAK' if okO else 'NIE'}")
V = {}
for v in ('old', 'new'):
    t = open(os.path.join(D, 'out', f'validation_{v}.txt')).read()
    m = re.search(r'Enforced checks:\s*(\d+)/(\d+)', t); e = re.search(r'EXIT (\d+)', t)
    V[v] = (int(m.group(1)), int(m.group(2)), int(e.group(1))) if m and e else None
    pr(f'  walidacja {v}: {V[v]}')
okV = V['new'] is not None and V['old'] is not None and V['new'][2] == 0 and V['new'][0] == V['new'][1] and V['new'][:2] == V['old'][:2]
pr(f"H795v (walidacja new EXIT 0, wszystkie kontrole, liczba jak old): {'TAK' if okV else 'NIE'}")
# --- H795t (unit check)
okT = True; nT = 0
for ln in open(os.path.join(D, 'out', 't795.txt')):
    w = ln.split(); g = lambda k, o=1: float(w[w.index(k) + o])
    e_new = float(w[w.index('nowy:') + 2]); l_new = float(w[w.index('nowy:') + 4]); L = float(w[w.index('L') + 1])
    okT &= abs(e_new) <= 1e-9 and abs(l_new/L - 1) <= 1e-9; nT += 1
pr(f"H795t (t795: {nT} linii, wodor i Ps, dipole 0 / rownolegle / antyrownolegle: |E_orb/E_S - 1| <= 1e-9, L = nominalne do 1e-9): {'TAK' if okT else 'NIE'}")
pr('  Kepler (CREM_KEPLER_START=1) E_orb/E_S - 1, bez dipoli:')
for ln in open(os.path.join(D, 'out', 't795.txt')):
    if ' dip 0 ' in ln: w = ln.split(); pr(f'    {w[0]} n {w[4]} L {w[6]}: {w[w.index("Kepler:") + 2]} | r_new/r_Kep - 1 {w[w.index("r_new/r_Kep-1") + 1]}')
okS = True
pr('=== S: wodor, linia TORQUE794 nowego measurementState (E0 z energia dipol-dipol) ===')
for n, L in TORI:
    a = relstart(f'{O}/S/F1_n{n}L{lab(L)}.err'); b = relstart(f'{O}/S/F0_n{n}L{lab(L)}.err')
    okS &= abs(a['E0']/a['ES'] - 1) <= 1e-4 and abs(b['E0']/b['ES'] - 1) <= 1e-4
    pr(f'  n {n} L {L:3}: F = 1 {a["E0"]/a["ES"] - 1:+.3e} | F = 0 {b["E0"]/b["ES"] - 1:+.3e}')
pr(f"H795t (S, wodor F = 0 i F = 1: |E0/E_S - 1| <= 1e-4): {'TAK' if okS else 'NIE'}")
okT0 = True
for n, L in TORI:
    a = relstart(f'{O}/S0/F1_n{n}L{lab(L)}.err'); b = relstart(f'{R794}/B/proton_n{n}L{lab(L)}_5.err')
    q = (a['E0']/a['ES'] - 1)/(b['E0']/b['ES'] - 1) - 1; okT0 &= abs(q) <= 1e-6
    pr(f'  S0 n {n} L {L:3}: {a["E0"]/a["ES"] - 1:+.4e} (794 {b["E0"]/b["ES"] - 1:+.4e}) -> {q:+.1e}')
pr(f"H795t0 (S0 = 794 do 1e-6 wzglednie): {'TAK' if okT0 else 'NIE'}")
# --- rule
okD = okM = True; out = {}; ALL = {}
pr('=== G (start produkcyjny, F = 1, CREM_NO_DIPOLE_FORCE) ===')
for t in ('proton', 'antiproton'):
    for s in '456':
        W = {k: omega(f'{O}/G/{t}_n{k[0]}L{lab(k[1])}_{s}.err') for k in TORI}
        W4 = {k: omega(f'{R794}/B/{t}_n{k[0]}L{lab(k[1])}_{s}.err') for k in TORI}
        for k in TORI: ALL.setdefault(k, []).append(W[k])
        K = {k: W[k]*k[0]**3*k[1]**2 for k in TORI}; sp = max(K.values())/min(K.values()) - 1
        pr(f'  {t} ziarno {s}: rozrzut Omega n^3 L^2 {sp*100:.3f} % | wobec 794 B: ' + ', '.join(f'n{k[0]}L{lab(k[1])} {W[k]/W4[k] - 1:+.1e}' for k in TORI))
        for n, l, st in ((2, 1, '2p'), (3, 1, '3p'), (3, 2, '3d')):
            F = lambda L: L*W[(n, L)]; od = F(l) - F(l + 1); q = od/W[(n, l + 0.5)]; qt = (l + 0.5)**2/(l*(l + 1)); okD &= abs(q/qt - 1) <= 0.01
            out.setdefault(st, []).append(q/qt)
            if st == '2p': okM &= abs(od/FS - 1) <= 0.01; out.setdefault('2P', []).append(od)
            pr(f'    {st}: {units(od)} | {q/qt:.5f} x regula' + (f' | {od/FS:.5f} pomiaru' if st == '2p' else ''))
pr(f"H795d (F = 1: 2p {min(out['2p']):.5f}-{max(out['2p']):.5f}, 3p {min(out['3p']):.5f}-{max(out['3p']):.5f}, 3d {min(out['3d']):.5f}-{max(out['3d']):.5f} w 1 %): {'TAK' if okD else 'NIE'}")
pr(f"H795m (F = 1: 2P {min(out['2P'])/1e6:.2f}-{max(out['2P'])/1e6:.2f} MHz w 1 % pomiaru): {'TAK' if okM else 'NIE'}")
M = {k: np.mean(v) for k, v in ALL.items()}; od = M[(2, 1)] - 2*M[(2, 2)]
pr(f'  srednie 6 biegow: 2P {units(od)} -> {od/FS:.5f} pomiaru (10 969,045(15) MHz, Hagley i Pipkin 1994)')
okF = True; mx = 0.0
for t in ('proton', 'antiproton'):
    for k in TORI:
        q = omega(f'{O}/G0/{t}_n{k[0]}L{lab(k[1])}_5.err')/omega(f'{O}/G/{t}_n{k[0]}L{lab(k[1])}_5.err') - 1; mx = max(mx, abs(q)); okF &= abs(q) <= 1e-4
pr(f"H795f (F = 0 / F = 1 ziarno 5, max {mx:.1e} <= 1e-4; tozsamosc, nie test singletu): {'TAK' if okF else 'NIE'}")
R21 = [x*9/8 for x in out['2p']]
pr(f'  o-Ps: 8/9 x R x 0,9998..0,9999 = {8/9*min(R21)*0.9998:.5f}-{8/9*max(R21)*0.9999:.5f} multipletu (Hagena i in., PRL 71, 2887 (1993)); p-Ps bez czesci spinowej')

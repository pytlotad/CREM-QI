# Audit 830: spin arm (charge-moment forces on), retarded (O) vs conservative (Z), apsidal rate n = 2, l = 1.  Usage: python3 ana830.py
# -> out/ana830.txt.  Written before the runs; gauss, read and nu_of copied from ana828.py (read gains half = first half of apocentres).
# Node u5c at 2048 taken from ../section-828/runs/runs828.tar.xz (D_* members).
import os, math, tarfile
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(D, 'out', 'runs'); C8 = os.path.join(D, 'out', 'c828')
HEV, C = 4.135667696e-15, 299792458.0
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:+.4f} MHz = {f*HEV*1e6:+.5f} ueV' + (f' (|lambda| {C/abs(f)*100:.4g} cm)' if f else '')
def gauss(N):
    x, w = np.polynomial.legendre.leggauss(N); u = 0.75 + 0.25*x[::-1]; return 1/u, 0.25*w[::-1]
NODES = dict(zip([f'u5{c}' for c in 'abcde'], zip(*gauss(5))))
def read(path, tag='mechB0', half=False):
    rows = []; has_sum = False
    for ln in open(path):
        if ln.startswith(f'PAIR TORQUE799 {tag} sum'): has_sum = True
        if not ln.startswith(f'PAIR APSE808 {tag} '): continue
        w = ln.split(); tc = float(w[w.index('tc') + 1])
        r = [np.array([float(v) for v in w[w.index(f'r{q}') + 1:w.index(f'r{q}') + 4]]) for q in range(3)]
        rows.append((tc, r))
    if half: rows = rows[:len(rows)//2]
    r0 = rows[0][1]; Lh = np.cross(r0[0], r0[2]); Lh /= np.linalg.norm(Lh)
    e1 = r0[1]/np.linalg.norm(r0[1]); e1 -= Lh*np.dot(Lh, e1); e1 /= np.linalg.norm(e1); e2 = np.cross(Lh, e1)
    T = []; P = []; nbad = 0; prev = None
    for tc, r in rows:   # = ana808.apse
        th = np.unwrap([math.atan2(np.dot(x, e2), np.dot(x, e1)) for x in r]); rho = [np.linalg.norm(x) for x in r]
        x = th - th[1]; Am = np.vstack([np.ones(3), x, x**2]).T; a0, b, cq = np.linalg.solve(Am, rho)
        xs = -b/(2*cq); hh = max(abs(x[0]), abs(x[2]))
        if not (cq < 0 and x[0] - 0.5*hh <= xs <= x[2] + 0.5*hh): nbad += 1
        phi = th[1] + xs
        if prev is not None: phi += 2*math.pi*round((prev - phi)/(2*math.pi))
        prev = phi; T.append(tc); P.append(phi)
    T = np.array(T); P = np.array(P)
    lin = np.polyfit(T, P, 1)[0]; q = np.polyfit(T, P, 2); om0 = q[1]
    Tk = np.diff(T); tm = 0.5*(T[1:] + T[:-1]); slope = np.polyfit(tm, np.log(Tk), 1)[0]
    return dict(om0=om0, lin=lin, slope=slope, n=len(T), nbad=nbad, has_sum=has_sum)
nu_of = lambda Om: sum(w*L**2*Om[k] for k, (L, w) in NODES.items())/(2*math.pi)
if not os.path.isdir(C8):
    with tarfile.open(os.path.join(D, '..', 'section-828', 'runs', 'runs828.tar.xz')) as tf:
        tf.extractall(C8, members=[m for m in tf.getmembers() if '/D_' in m.name and 'u5c_s2048' in m.name])
SYS = {'ops': ('F1_ops', 'o-Ps (PHEN=2)'), 'pps': ('F0_pps', 'p-Ps (PHEN=1)'), 'H1': ('F1_H', 'wodor F = 1 (PHEN=2)'), 'H0': ('F0_H', 'wodor F = 0 (PHEN=1)')}
CH828 = {'ops': -0.1557e6, 'pps': -0.1557e6, 'H1': 0.0001e6, 'H0': 0.0001e6}   # 828 charge arm Delta (2048), MHz -> Hz
path = lambda m, nm, k, h=2048: (f'{C8}/runs/{m}/D_{nm}_{k}_s{h}.err' if (k == 'u5c' and h == 2048) else f'{R}/{m}/D_{nm}_{k}_s{h}.err')
V = {}; okF = True
for m in ('Z', 'O'):
    for sk, (nm, _) in SYS.items():
        for k in NODES:
            for hf in (False, True):
                v = read(path(m, nm, k), half=hf); V[(m, sk, k, hf)] = v
                if not hf: okF &= v['has_sum'] and v['n'] >= 195
    for sk in ('ops', 'H0'):
        V[(m, sk, 'u5c', 4096)] = v = read(path(m, SYS[sk][0], 'u5c', 4096)); okF &= v['has_sum'] and v['n'] >= 195
pr(f"  pliki: {sum(1 for k in V if k[3] is not True)}; wszystkie z linia sum i >= 195 apocentrami: {okF}; zlych wierzcholkow {sum(v['nbad'] for k, v in V.items() if k[3] is not True)}")
same = open(f'{R}/Z/D_F1_ops_u5c_s2048.err', 'rb').read() == open(path('Z', 'F1_ops', 'u5c'), 'rb').read()
pr(f"H830b (powtorzony Z/D_F1_ops_u5c_s2048 = 828 bajt w bajt): {'TAK' if same else 'NIE'}")
Dl = {}; Dh = {}
for sk, (nm, lab) in SYS.items():
    nz = nu_of({k: V[('Z', sk, k, False)]['om0'] for k in NODES}); no = nu_of({k: V[('O', sk, k, False)]['om0'] for k in NODES})
    hz = nu_of({k: V[('Z', sk, k, True)]['om0'] for k in NODES}); ho = nu_of({k: V[('O', sk, k, True)]['om0'] for k in NODES})
    Dl[sk] = no - nz; Dh[sk] = ho - hz
    pr(f'  {lab}: nu_orb Z {nz/1e6:.4f} MHz, O {no/1e6:.4f} MHz; Delta = O - Z = {units(no - nz)}; wzglednie {(no - nz)/nz:+.2e}')
    pr(f'    polowa okna: Delta {units(ho - hz)}; czesc spinowa Delta - Delta_ladunkowe(828) {units(no - nz - CH828[sk])}')
    pr('    wezly dOmega_0/Omega_0: ' + ', '.join(f"{k} {V[('O', sk, k, False)]['om0']/V[('Z', sk, k, False)]['om0'] - 1:+.2e}" for k in NODES))
    pr('    nachylenie ln T  Z: ' + ', '.join(f"{V[('Z', sk, k, False)]['slope']:+.2e}" for k in NODES) + ';  O: ' + ', '.join(f"{V[('O', sk, k, False)]['slope']:+.2e}" for k in NODES))
okK = True
for sk in ('ops', 'H0'):
    a = V[('O', sk, 'u5c', False)]['om0']/V[('Z', sk, 'u5c', False)]['om0'] - 1; b = V[('O', sk, 'u5c', 4096)]['om0']/V[('Z', sk, 'u5c', 4096)]['om0'] - 1
    okK &= abs(b - a) <= 0.1*abs(a); pr(f'  krok u5c {SYS[sk][1]}: dOmega_0/Omega_0 2048 {a:+.3e}, 4096 {b:+.3e}')
pr(f"H830k (krok: |4096 - 2048| <= 0,1 |2048|, o-Ps i wodor F = 0): {'TAK' if okK else 'NIE'}")
okW = all(abs(Dh[sk] - Dl[sk]) <= max(0.2*abs(Dl[sk]), 0.05e6) for sk in SYS)
pr(f"H830w (okno: |Delta(polowa) - Delta| <= max(0,2 |Delta|, 0,05 MHz), cztery uklady): {'TAK' if okW else 'NIE'}")
pr(f"H830o (GLOWNA: o-Ps |Delta| >= 0,61 MHz): {'TAK' if abs(Dl['ops']) >= 0.61e6 else 'NIE'}")
pr(f"H830s (p-Ps |Delta| < 0,61 MHz): {'TAK' if abs(Dl['pps']) < 0.61e6 else 'NIE'}")
pr(f"H830h (wodor F = 0 i F = 1: |Delta| >= 0,015 MHz): {'TAK' if abs(Dl['H0']) >= 0.015e6 and abs(Dl['H1']) >= 0.015e6 else 'NIE'}")

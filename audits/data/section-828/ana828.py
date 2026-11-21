# Audit 828: retarded (O) vs conservative Coulomb + Darwin (Z) apsidal rate at n = 2, l = 1.  Usage: python3 ana828.py -> out/ana828.txt.
# Written before the runs; read() and nu_of() as ana813.py.
import os, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(D, 'out', 'runs')
HEV, C = 4.135667696e-15, 299792458.0
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:+.4f} MHz = {f*HEV*1e6:+.5f} ueV' + (f' (|lambda| {C/abs(f)*100:.4g} cm)' if f else '')
def gauss(N):
    x, w = np.polynomial.legendre.leggauss(N); u = 0.75 + 0.25*x[::-1]; return 1/u, 0.25*w[::-1]
NODES = dict(zip([f'u5{c}' for c in 'abcde'], zip(*gauss(5))))
def read(path, tag='mechB0'):
    rows = []; has_sum = False
    for ln in open(path):
        if ln.startswith(f'PAIR TORQUE799 {tag} sum'): has_sum = True
        if not ln.startswith(f'PAIR APSE808 {tag} '): continue
        w = ln.split(); tc = float(w[w.index('tc') + 1])
        r = [np.array([float(v) for v in w[w.index(f'r{q}') + 1:w.index(f'r{q}') + 4]]) for q in range(3)]
        rows.append((tc, r))
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
V = {}; okF = True
for m in ('Z', 'O'):
    for k in NODES:
        for key, fn in ((('ops', 2048), f'F1_ops_{k}_s2048'), (('ops', 4096), f'F1_ops_{k}_s4096'), (('pps', 2048), f'F0_pps_{k}_s2048'),
                        (('H1', 2048), f'F1_H_{k}_s2048'), (('H0', 2048), f'F0_H_{k}_s2048')):
            v = read(f'{R}/{m}/{fn}.err'); V[(m, key, k)] = v; okF &= v['has_sum'] and v['n'] >= 195
    for sk in ('ops', 'pps', 'H1', 'H0'):
        nm = {'ops': 'F1_ops', 'pps': 'F0_pps', 'H1': 'F1_H', 'H0': 'F0_H'}[sk]
        v = read(f'{R}/{m}/D_{nm}_u5c_s2048.err'); V[(m, ('D' + sk, 2048), 'u5c')] = v; okF &= v['has_sum'] and v['n'] >= 195
pr(f'  pliki: {len(V)}; wszystkie z linia sum i >= 195 apocentrami: {okF}; zlych wierzcholkow {sum(v["nbad"] for v in V.values())}')
LAB = {('ops', 2048): 'o-Ps (PHEN=2) 2048', ('ops', 4096): 'o-Ps (PHEN=2) 4096', ('pps', 2048): 'p-Ps (PHEN=1) 2048',
       ('H1', 2048): 'wodor F = 1 (PHEN=2) 2048', ('H0', 2048): 'wodor F = 0 (PHEN=1) 2048'}
Dl = {}
for key, lab in LAB.items():
    nz = nu_of({k: V[('Z', key, k)]['om0'] for k in NODES}); no = nu_of({k: V[('O', key, k)]['om0'] for k in NODES}); Dl[key] = no - nz
    pr(f'  {lab}: nu_orb Z {nz/1e6:.4f} MHz, O {no/1e6:.4f} MHz; Delta = O - Z = {units(no - nz)}; wzglednie {(no - nz)/nz:+.2e}')
    pr('    wezly dOmega_0/Omega_0: ' + ', '.join(f"{k} {V[('O', key, k)]['om0']/V[('Z', key, k)]['om0'] - 1:+.2e}" for k in NODES)
       + '; nachylenie ln T O: ' + ', '.join(f"{V[('O', key, k)]['slope']:+.2e}" for k in NODES))
dk = abs(Dl[('ops', 4096)] - Dl[('ops', 2048)])
pr(f"H828k (|Delta(4096) - Delta(2048)| o-Ps <= 0,010 MHz): {'TAK' if dk <= 0.010e6 else 'NIE'} ({dk/1e6:.4f} MHz)")
dc = max(abs(V[(m, ('pps', 2048), k)]['om0']/V[(m, ('ops', 2048), k)]['om0'] - 1) for m in ('Z', 'O') for k in NODES)
pr(f"H828c (o-Ps = p-Ps do 1e-8, ramie ladunkowe): {'TAK' if dc <= 1e-8 else 'NIE'} (max {dc:.2e})")
okP = abs(Dl[('ops', 4096)]) < 0.61e6 and abs(Dl[('ops', 2048)]) < 0.61e6 and abs(Dl[('pps', 2048)]) < 0.61e6
pr(f"H828p (GLOWNA: |Delta| < 0,61 MHz, o-Ps i p-Ps): {'TAK' if okP else 'NIE'}")
okH = abs(Dl[('H1', 2048)]) < 0.015e6 and abs(Dl[('H0', 2048)]) < 0.015e6
pr(f"H828h (wodor F = 0 i F = 1: |Delta| < 0,015 MHz): {'TAK' if okH else 'NIE'}")
pr('  RAPORT ramie spinowe (u5c, 2048, sily ladunek-moment wlaczone): dOmega_0/Omega_0 O - Z; obok ramie ladunkowe u5c')
for sk, lab in (('ops', 'o-Ps (PHEN=2)'), ('pps', 'p-Ps (PHEN=1)'), ('H1', 'wodor F = 1 (PHEN=2)'), ('H0', 'wodor F = 0 (PHEN=1)')):
    a = V[('O', ('D' + sk, 2048), 'u5c')]['om0']/V[('Z', ('D' + sk, 2048), 'u5c')]['om0'] - 1
    b = V[('O', (sk, 2048), 'u5c')]['om0']/V[('Z', (sk, 2048), 'u5c')]['om0'] - 1
    pr(f'    {lab}: spinowe {a:+.3e}, ladunkowe {b:+.3e}, roznica {a - b:+.3e}; Omega_0 Z spin / bez spinu - 1 '
       f"{V[('Z', ('D' + sk, 2048), 'u5c')]['om0']/V[('Z', (sk, 2048), 'u5c')]['om0'] - 1:+.3e}")

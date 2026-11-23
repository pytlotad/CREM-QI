# Audit 836: retarded mode O started at the Tulczyjew centroid (OS) vs Z and O (830/828).  Usage: python3 ana836.py -> out/ana836.txt.
# Written before the runs; gauss, read, nu_of copied from ana830.py.
import os, math, tarfile
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(D, 'out', 'runs'); BB = os.path.join(D, 'out', 'base')
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
    for tc, r in rows:
        th = np.unwrap([math.atan2(np.dot(x, e2), np.dot(x, e1)) for x in r]); rho = [np.linalg.norm(x) for x in r]
        x = th - th[1]; Am = np.vstack([np.ones(3), x, x**2]).T; a0, b, cq = np.linalg.solve(Am, rho)
        xs = -b/(2*cq); hh = max(abs(x[0]), abs(x[2]))
        if not (cq < 0 and x[0] - 0.5*hh <= xs <= x[2] + 0.5*hh): nbad += 1
        phi = th[1] + xs
        if prev is not None: phi += 2*math.pi*round((prev - phi)/(2*math.pi))
        prev = phi; T.append(tc); P.append(phi)
    T = np.array(T); P = np.array(P)
    q = np.polyfit(T, P, 2); om0 = q[1]
    Tk = np.diff(T); tm = 0.5*(T[1:] + T[:-1]); slope = np.polyfit(tm, np.log(Tk), 1)[0]
    return dict(om0=om0, slope=slope, n=len(T), nbad=nbad, has_sum=has_sum)
nu_of = lambda Om: sum(w*L**2*Om[k] for k, (L, w) in NODES.items())/(2*math.pi)
for sec, filt in (('830', lambda n: '/D_' in n), ('828', lambda n: '/D_' in n and 'u5c_s2048' in n)):
    dst = os.path.join(BB, sec)
    if not os.path.isdir(dst):
        with tarfile.open(os.path.join(D, '..', f'section-{sec}', 'runs', f'runs{sec}.tar.xz')) as tf:
            tf.extractall(dst, members=[m for m in tf.getmembers() if filt(m.name)])
base = lambda m, nm, k: f"{BB}/{'828' if k == 'u5c' else '830'}/runs/{m}/D_{nm}_{k}_s2048.err"
SYS = {'ops': ('F1_ops', 'o-Ps (PHEN=2)'), 'H0': ('F0_H', 'wodor F = 0 (PHEN=1)'), 'H1': ('F1_H', 'wodor F = 1 (PHEN=2)'), 'pps': ('F0_pps', 'p-Ps (PHEN=1)')}
DZ = {'ops': 3.1091e6, 'H0': 9.9695e6, 'H1': -2.0140e6, 'pps': 0.0026e6}    # delta_Z (834)
DS = {'ops': 0.0411e6, 'H0': 0.3404e6, 'H1': -0.2155e6, 'pps': 0.0026e6}    # delta_S (835)
CH = {'ops': -0.1557e6, 'H0': 0.0001e6, 'H1': 0.0001e6, 'pps': -0.1557e6}  # charge arm (828)
same = open(f'{R}/K/D_F1_ops_u5c_s2048.err', 'rb').read() == open(base('O', 'F1_ops', 'u5c'), 'rb').read()
pr(f"H836b (K = O/D_F1_ops_u5c_s2048 828 bajt w bajt): {'TAK' if same else 'NIE'}")
okM = True; okF = True
for sk, (nm, lab) in SYS.items():
    V = {}
    for m in ('Z', 'O', 'OS'):
        for hf in (False, True):
            V[(m, hf)] = {k: read(f'{R}/OS/D_{nm}_{k}_s2048.err' if m == 'OS' else base(m, nm, k), half=hf) for k in NODES}
    okF &= all(v['has_sum'] and v['n'] >= 195 for v in V[('OS', False)].values())
    nu = {key: nu_of({k: v['om0'] for k, v in d.items()}) for key, d in V.items()}
    d830 = nu[('O', False)] - nu[('Z', False)]; dS = nu[('OS', False)] - nu[('Z', False)]; dSh = nu[('OS', True)] - nu[('Z', True)]
    rem = d830 - dS; exp = DZ[sk] - DS[sk]
    pr(f'  {lab}: Delta_830 = O - Z {units(d830)}')
    pr(f'    Delta_S = OS - Z {units(dS)};  polowa okna {units(dSh)}')
    pr(f'    zabrane przesunieciem startu Delta_830 - Delta_S {units(rem)}; oczekiwane (delta_Z - delta_S, 834/835) {units(exp)}'
       + (f'; |roznica| / |delta_Z| {abs(rem - exp)/abs(DZ[sk]):.3f}' if sk != 'pps' else ''))
    pr(f'    Delta_S - ladunkowe(828) {units(dS - CH[sk])}')
    pr('    wezly dOmega/Omega (OS - Z): ' + ', '.join(f"{k} {V[('OS', False)][k]['om0']/V[('Z', False)][k]['om0'] - 1:+.3e}" for k in NODES)
       + f"; zle apsydy OS {[V[('OS', False)][k]['nbad'] for k in NODES]}")
    pr('    nachylenie ln T  O: ' + ', '.join(f"{V[('O', False)][k]['slope']:+.2e}" for k in NODES) + ';  OS: ' + ', '.join(f"{V[('OS', False)][k]['slope']:+.2e}" for k in NODES))
    if sk != 'pps': okM &= abs(rem - exp) <= 0.15*abs(DZ[sk])
pr(f'  pliki OS: wszystkie z linia sum i >= 195 apocentrami: {okF}')
pr(f"H836m (GLOWNA: |(Delta_830 - Delta_S) - (delta_Z - delta_S)| <= 0,15 |delta_Z|, o-Ps, F = 0, F = 1): {'TAK' if okM else 'NIE'}")

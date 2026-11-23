# Audit 832: retarded spin-arm parts (Mq, qM, MM, T) replaced one at a time by their instantaneous forms.  Usage: python3 ana832.py
# -> out/ana832.txt.  Written before the runs; gauss, read, nu_of copied from ana830.py.  O and Z from ../section-830 (u5a, u5b, u5d, u5e)
# and ../section-828 (u5c) archives, extracted into out/base.
import os, math, tarfile
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(D, 'out', 'runs'); BS = os.path.join(D, 'out', 'base')
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
if not os.path.isdir(BS):
    for sec in ('830', '828'):
        with tarfile.open(os.path.join(D, '..', f'section-{sec}', 'runs', f'runs{sec}.tar.xz')) as tf:
            ms = [m for m in tf.getmembers() if '/D_' in m.name and m.name.endswith('_s2048.err') and (('u5c' in m.name) == (sec == '828'))]
            tf.extractall(os.path.join(BS, sec), members=ms)
base = lambda m, nm, k: os.path.join(BS, '828' if k == 'u5c' else '830', 'runs', m, f'D_{nm}_{k}_s2048.err')
SYS = {'ops': ('F1_ops', 'o-Ps (PHEN=2)'), 'H0': ('F0_H', 'wodor F = 0 (PHEN=1)'), 'H1': ('F1_H', 'wodor F = 1 (PHEN=2)')}
DCH = {'ops': -0.1557e6, 'H0': 0.0001e6, 'H1': 0.0001e6}          # 828 charge arm Delta, Hz
RCH = {'ops': -2.549e-5, 'pps': -2.549e-5, 'H0': 4.624e-8, 'H1': 4.612e-8}   # 828 charge arm u5c dOmega/Omega
PARTS = ('Mq', 'qM', 'MM', 'T')
same = open(f'{R}/K/D_F1_ops_u5c_s2048.err', 'rb').read() == open(base('O', 'F1_ops', 'u5c'), 'rb').read()
pr(f"H832b (K = 828 bajt w bajt): {'TAK' if same else 'NIE'}")
okF = True; okA = True; frac = {}
for sk, (nm, lab) in SYS.items():
    V = {m: {k: read(base(m, nm, k)) for k in NODES} for m in ('O', 'Z')}
    for p in PARTS:
        V[p] = {k: read(f'{R}/{p}/D_{nm}_{k}_s2048.err') for k in NODES}
    for m in V:
        for k in NODES: okF &= V[m][k]['has_sum'] and V[m][k]['n'] >= 195
    nu = {m: nu_of({k: V[m][k]['om0'] for k in NODES}) for m in V}
    Ds = nu['O'] - nu['Z'] - DCH[sk]
    pr(f'  {lab}: Delta = O - Z {units(nu["O"] - nu["Z"])}; czesc spinowa Delta_s {units(Ds)}')
    tot = 0.0
    for p in PARTS:
        d = nu['O'] - nu[p]; tot += d; frac[(sk, p)] = d/Ds
        pr(f'    delta_{p} = O - P_{p}: {units(d)}; ulamek Delta_s {d/Ds:+.3f}; wezly dOmega/Omega '
           + ', '.join(f"{k} {V['O'][k]['om0']/V[p][k]['om0'] - 1:+.2e}" for k in NODES))
    pr(f'    suma delta {units(tot)}; suma / Delta_s {tot/Ds:+.3f}')
    okA &= abs(tot - Ds) <= 0.2*abs(Ds)
pr(f'  pliki: wszystkie z linia sum i >= 195 apocentrami: {okF}')
pr(f"H832a (|sum delta - Delta_s| <= 0,2 |Delta_s|, o-Ps, wodor F = 0, F = 1): {'TAK' if okA else 'NIE'}")
okC = True
for sk, nm, lab in (('ops', 'F1_ops', 'o-Ps (PHEN=2)'), ('pps', 'F0_pps', 'p-Ps (PHEN=1)'), ('H0', 'F0_H', 'wodor F = 0 (PHEN=1)'), ('H1', 'F1_H', 'wodor F = 1 (PHEN=2)')):
    z = read(base('Z', nm, 'u5c'))['om0']; o = read(base('O', nm, 'u5c'))['om0']; a = read(f'{R}/A/D_{nm}_u5c_s2048.err')['om0']
    rA = a/z - 1; rO = o/z - 1; ok = abs(rA - RCH[sk]) <= 0.1*abs(rO - RCH[sk])
    if sk != 'pps': okC &= ok
    pr(f'  u5c {lab}: r_A {rA:+.3e}, r_lad (828) {RCH[sk]:+.3e}, r_O {rO:+.3e}; |r_A - r_lad| / |r_O - r_lad| {abs(rA - RCH[sk])/abs(rO - RCH[sk]):.3f}')
pr(f"H832c (domkniecie u5c: |r_A - r_lad| <= 0,1 |r_O - r_lad|, o-Ps, wodor F = 0, F = 1): {'TAK' if okC else 'NIE'}")
dom = {sk: max(PARTS, key=lambda p: frac[(sk, p)]) for sk in SYS}
okD = len(set(dom.values())) == 1 and all(frac[(sk, dom[sk])] >= 0.5 for sk in SYS)
pr('  czlon dominujacy: ' + ', '.join(f'{SYS[sk][1]} {dom[sk]} ({frac[(sk, dom[sk])]:+.3f})' for sk in SYS))
pr(f"H832d (GLOWNA: jeden czlon >= 50 % Delta_s, ten sam w o-Ps, wodorze F = 0 i F = 1): {'TAK' if okD else 'NIE'}")

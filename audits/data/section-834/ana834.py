# Audit 834: Tulczyjew Thomas momentum in the instantaneous sum (mode Z).  Usage: python3 ana834.py -> out/ana834.txt.  Written before
# the runs; gauss, read, nu_of, base extraction copied from ana832.py.  Z from ../section-830 (u5a, b, d, e) and ../section-828 (u5c).
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
SYS = {'ops': ('F1_ops', 'o-Ps (PHEN=2)'), 'H0': ('F0_H', 'wodor F = 0 (PHEN=1)'), 'H1': ('F1_H', 'wodor F = 1 (PHEN=2)'), 'pps': ('F0_pps', 'p-Ps (PHEN=1)')}
DT = {'ops': 2.9540e6, 'H0': 9.6536e6, 'H1': -1.6983e6}   # 832 delta_T, Hz
same = open(f'{R}/K/D_F1_ops_u5c_s2048.err', 'rb').read() == open(base('Z', 'F1_ops', 'u5c'), 'rb').read()
pr(f"H834b (K = Z 828 bajt w bajt): {'TAK' if same else 'NIE'}")
dZ = {}; okF = True
for sk, (nm, lab) in SYS.items():
    V = {'Z': {k: read(base('Z', nm, k)) for k in NODES}, 'ZT': {k: read(f'{R}/ZT/D_{nm}_{k}_s2048.err') for k in NODES}}
    for m in V:
        for k in NODES: okF &= V[m][k]['has_sum'] and V[m][k]['n'] >= 195
    nu = {m: nu_of({k: V[m][k]['om0'] for k in NODES}) for m in V}
    dZ[sk] = nu['ZT'] - nu['Z']
    pr(f'  {lab}: delta_Z = ZT - Z {units(dZ[sk])}' + (f'; delta_T (832) {units(DT[sk])}; delta_Z / delta_T {dZ[sk]/DT[sk]:+.3f}' if sk in DT else ''))
    pr('    wezly dOmega/Omega ' + ', '.join(f"{k} {(V['ZT'][k]['om0'] - V['Z'][k]['om0'])/V['Z'][k]['om0']:+.3e}" for k in NODES)
       + '; zle apsydy ZT ' + str([V['ZT'][k]['nbad'] for k in NODES]))
pr(f'  pliki: wszystkie z linia sum i >= 195 apocentrami: {okF}')
M3 = ('ops', 'H0', 'H1')
form = all(abs(dZ[s] - DT[s]) <= 0.25*abs(DT[s]) for s in M3); ret = all(abs(dZ[s]) <= 0.25*abs(DT[s]) for s in M3)
pr(f"H834d (GLOWNA: |delta_Z - delta_T| <= 0,25 |delta_T|): {'TAK' if form else 'NIE'}; opoznienie (|delta_Z| <= 0,25 |delta_T|): {'TAK' if ret else 'NIE'}"
   + ('' if form or ret else ' -- mieszane'))
ev = (dZ['H0'] + dZ['H1'])/2; ev832 = (DT['H0'] + DT['H1'])/2
pr(f"H834w (czesc parzysta wodoru {units(ev)} wobec 832 {units(ev832)}, do 25 %): {'TAK' if abs(ev - ev832) <= 0.25*abs(ev832) else 'NIE'}")
pr(f"  czesc nieparzysta wodoru (F0 - F1)/2: {units((dZ['H0'] - dZ['H1'])/2)} wobec 832 {units((DT['H0'] - DT['H1'])/2)}")

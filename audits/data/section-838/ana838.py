# Audit 838: start at the apocentre vs periapsis -- is the 835 residual start-independent?  Usage: python3 ana838.py -> out/ana838.txt.
# Written before the runs; gauss, read, nu_of copied from ana835.py.  Periapsis values from ana835 (out/ana835.txt numbers).
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
if not os.path.isdir(os.path.join(BS, '835')):
    with tarfile.open(os.path.join(D, '..', 'section-835', 'runs', 'runs835.tar.xz')) as tf:
        tf.extractall(os.path.join(BS, '835'), members=[m for m in tf.getmembers() if m.name.endswith('ZS/D_F1_ops_u5c_s2048.err')])
same = open(f'{R}/K/D_F1_ops_u5c_s2048.err', 'rb').read() == open(os.path.join(BS, '835', 'runs', 'ZS', 'D_F1_ops_u5c_s2048.err'), 'rb').read()
pr(f"H838b (K = ZS 835 bajt w bajt): {'TAK' if same else 'NIE'}")
for ln in open(f'{R}/ZA/D_F0_H_u5a_s2048.err'):
    if 'START838' in ln: pr('  ' + ln.strip())
SYS = {'ops': ('F1_ops', 'o-Ps (PHEN=2)'), 'H0': ('F0_H', 'wodor F = 0 (PHEN=1)'), 'H1': ('F1_H', 'wodor F = 1 (PHEN=2)')}
PZ = {'ops': 3.1091e6, 'H0': 9.9695e6, 'H1': -2.0140e6}; PSr = {'ops': 0.0411e6, 'H0': 0.3404e6, 'H1': -0.2155e6}
okK = okD = okF = True; dza = {}; dsa = {}
for sk, (nm, lab) in SYS.items():
    V = {m: {k: read(f'{R}/{m}/D_{nm}_{k}_s2048.err') for k in NODES} for m in ('ZA', 'ZTA', 'ZSA')}
    okF &= all(v['has_sum'] and v['n'] >= 195 and v['nbad'] == 0 for d in V.values() for v in d.values())
    nu = {m: nu_of({k: V[m][k]['om0'] for k in NODES}) for m in V}
    dza[sk] = nu['ZTA'] - nu['ZA']; dsa[sk] = nu['ZSA'] - nu['ZA']
    okK &= abs(dza[sk]) <= 0.5*abs(PZ[sk]); okD &= abs(dsa[sk] - PSr[sk]) <= 0.25*abs(PSr[sk])
    pr(f'  {lab}: nu(ZA) {nu["ZA"]/1e6:.4f} MHz')
    pr(f'    delta_Z(apo) = ZTA - ZA {units(dza[sk])}  (peri {PZ[sk]/1e6:+.4f}; apo / peri {dza[sk]/PZ[sk]:+.4f})')
    pr(f'    delta_S(apo) = ZSA - ZA {units(dsa[sk])}  (peri {PSr[sk]/1e6:+.4f}; apo / peri {dsa[sk]/PSr[sk]:+.4f})')
    pr(f'    |delta_Z(apo) - delta_S(peri)| / |delta_S(peri)| {abs(dza[sk] - PSr[sk])/abs(PSr[sk]):.3f}')
    for m in ('ZTA', 'ZSA'):
        pr(f'    wezly dOmega/Omega ({m} - ZA): ' + ', '.join(f"{k} {V[m][k]['om0']/V['ZA'][k]['om0'] - 1:+.3e}" for k in NODES))
    pr('    apocentra ZA/ZTA/ZSA: ' + ', '.join(f"{k} {V['ZA'][k]['n']}/{V['ZTA'][k]['n']}/{V['ZSA'][k]['n']}" for k in NODES))
pr(f'  pliki: wszystkie z linia sum, >= 195 apocentrami i 0 zlych apsyd: {okF}')
pr(f"H838k (|delta_Z(apo)| <= 0,5 |delta_Z(peri)|): {'TAK' if okK else 'NIE'}")
pr(f"H838d (GLOWNA: |delta_S(apo) - delta_S(peri)| <= 0,25 |delta_S(peri)|): {'TAK' if okD else 'NIE'}")
for nmq, d in (('delta_Z(apo)', dza), ('delta_S(apo)', dsa)):
    pr(f'  wodor {nmq}: czesc parzysta {(d["H0"] + d["H1"])/2/1e6:+.4f}, nieparzysta (F1 - F0)/2 {(d["H1"] - d["H0"])/2/1e6:+.4f} MHz')

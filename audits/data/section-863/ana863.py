# Audit 863: shift of the default p814 measurement (ZDC -> O) from archives: retarded OS 836 - OZ 846, conservative ZS 835 - ZDC 846, R = O - Z.
# Written after the controls (report prescribed in prereg863); header (units, gauss, read, nu_of) copied from ana862.py.
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
BB = os.path.join(D, 'out', 'base'); os.makedirs(BB, exist_ok=True)
for sec in ('835', '836', '846'):
    with tarfile.open(os.path.join(D, '..', f'section-{sec}', 'runs', f'runs{sec}.tar.xz')) as t: t.extractall(os.path.join(BB, sec), filter='data')
P = {'OS': '836/runs/OS', 'ZS': '835/runs/ZS', 'OZ': '846/runs/OZ', 'ZDC': '846/runs/ZDC'}
CH = (('o-Ps (PHEN=2)', 'D_F1_ops'), ('p-Ps (PHEN=1)', 'D_F0_pps'), ('wodor F = 0 (PHEN=1)', 'D_F0_H'), ('wodor F = 1 (PHEN=2)', 'D_F1_H'))
for lab, stem in CH:
    nu = {k: nu_of({n: read(f'{BB}/{p}/{stem}_{n}_s2048.err')['om0'] for n in NODES}) for k, p in P.items()}
    pr(f'{lab}:')
    for k in P: pr(f'  nu({k:3s}) = {units(nu[k])}')
    pr(f'  przesuniecie sumy opoznionej ZDC -> O   (OS - OZ)  {units(nu["OS"] - nu["OZ"])}')
    pr(f'  przesuniecie sumy natychmiastowej ZDC -> O (ZS - ZDC) {units(nu["ZS"] - nu["ZDC"])}')
    pr(f'  R domyslne nowe (OS - ZS) {units(nu["OS"] - nu["ZS"])};  R 848-862 (OZ - ZDC) {units(nu["OZ"] - nu["ZDC"])}')
import shutil; shutil.rmtree(BB)

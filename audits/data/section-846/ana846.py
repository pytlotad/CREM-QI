# Audit 846: retardation share in the ZDC form, R_ZDC = nu(OZ) - nu(ZDC), vs R_old = nu(OS 836) - nu(ZS 835).
# Written before the runs; header (gauss, read, nu_of) copied from ana836.py.
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

for sec, sub in (('835', '/ZS/'), ('836', '/OS/')):
    dst = os.path.join(BB, sec)
    if not os.path.isdir(dst):
        with tarfile.open(os.path.join(D, '..', f'section-{sec}', 'runs', f'runs{sec}.tar.xz')) as tf:
            tf.extractall(dst, members=[m for m in tf.getmembers() if sub in m.name])
SYS = {'ops': ('F1_ops', 'o-Ps (PHEN=2)'), 'pps': ('F0_pps', 'p-Ps (PHEN=1)'), 'H0': ('F0_H', 'wodor F = 0 (PHEN=1)'), 'H1': ('F1_H', 'wodor F = 1 (PHEN=2)')}
txt = lambda p: open(p).read()
same = open(f'{R}/K/D_F1_ops_u5c_s2048.err', 'rb').read() == open(f'{BB}/836/runs/OS/D_F1_ops_u5c_s2048.err', 'rb').read()
lines = all('PAIR ORETZDC846 on' in txt(f'{R}/OZ/D_{nm}_{k}_s2048.err') for nm, _ in SYS.values() for k in NODES) \
    and all('PAIR SODIRAC841 on' in txt(f'{R}/ZDC/D_{nm}_{k}_s2048.err') and 'PAIR SOCANON843 on' in txt(f'{R}/ZDC/D_{nm}_{k}_s2048.err') for nm, _ in SYS.values() for k in NODES) \
    and 'ORETZDC846' not in txt(f'{R}/K/D_F1_ops_u5c_s2048.err')
pr(f"H846b (K = OS 836 bajt w bajt; ORETZDC846 we wszystkich OZ, brak w K; SODIRAC841 + SOCANON843 we wszystkich ZDC): {'TAK' if same and lines else 'NIE'}  (K {same}, wiersze {lines})")
okM = True; RES = {}
for sk, (nm, lab) in SYS.items():
    V = {}
    for m, path in (('OZ', f'{R}/OZ'), ('ZDC', f'{R}/ZDC'), ('OS', f'{BB}/836/runs/OS'), ('ZS', f'{BB}/835/runs/ZS')):
        for hf in (False, True):
            V[(m, hf)] = {k: read(f'{path}/D_{nm}_{k}_s2048.err', half=hf) for k in NODES}
    nu = {key: nu_of({k: v['om0'] for k, v in d.items()}) for key, d in V.items()}
    Rz = nu[('OZ', False)] - nu[('ZDC', False)]; Ro = nu[('OS', False)] - nu[('ZS', False)]
    Rzh = nu[('OZ', True)] - nu[('ZDC', True)]; Roh = nu[('OS', True)] - nu[('ZS', True)]
    tol = max(0.05e6, 0.10*abs(Ro)); ok = abs(Rz - Ro) <= tol; okM &= ok; RES[sk] = (Rz, Ro)
    pr(f'  {lab}: R_ZDC = OZ - ZDC {units(Rz)}')
    pr(f'    R_old = OS - ZS {units(Ro)};  R_ZDC - R_old {(Rz - Ro)/1e6:+.4f} MHz (granica {tol/1e6:.4f}) -> {"TAK" if ok else "NIE"}')
    pr(f'    polowa okna: R_ZDC {Rzh/1e6:+.4f}, R_old {Roh/1e6:+.4f} MHz;  OZ - OS {(nu[("OZ", False)] - nu[("OS", False)])/1e6:+.4f}, ZDC - ZS {(nu[("ZDC", False)] - nu[("ZS", False)])/1e6:+.4f} MHz')
    pr('    wezly R_ZDC / R_old: ' + ', '.join(f"{k} {(V[('OZ', False)][k]['om0'] - V[('ZDC', False)][k]['om0'])/(V[('OS', False)][k]['om0'] - V[('ZS', False)][k]['om0']):.3f}" for k in NODES))
    pr(f"    apocentra min OZ {min(v['n'] for v in V[('OZ', False)].values())}, ZDC {min(v['n'] for v in V[('ZDC', False)].values())}; zle apsydy OZ {sum(v['nbad'] for v in V[('OZ', False)].values())}, ZDC {sum(v['nbad'] for v in V[('ZDC', False)].values())}")
    pr('    nachylenie ln T  OZ: ' + ', '.join(f"{V[('OZ', False)][k]['slope']:+.2e}" for k in NODES) + ';  OS: ' + ', '.join(f"{V[('OS', False)][k]['slope']:+.2e}" for k in NODES))
for tag, i in (('R_ZDC', 0), ('R_old', 1)):
    ev, od = (RES['H0'][i] + RES['H1'][i])/2, (RES['H1'][i] - RES['H0'][i])/2
    pr(f'  wodor {tag}: parzysta {units(ev)}, nieparzysta (F1 - F0)/2 {units(od)}')
pr(f"H846m (GLOWNA: |R_ZDC - R_old| <= max(0,05 MHz; 0,10 |R_old|) w o-Ps, p-Ps, F = 0, F = 1): {'TAK' if okM else 'NIE'}")

# Audit 862: R_ZDC with the Tulczyjew start, R_ZDC,T = nu(OZT) - nu(ZDCT), vs R_old (OS 836 - ZS 835) and R_ZDC (846).
# Written before the runs; header (units, gauss, read, nu_of) copied from ana846.py.
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


ARCH = {'846': ('/OZ/', '/ZDC/'), '835': ('/ZS/',), '836': ('/OS/',)}
for sec, subs in ARCH.items():
    dst = os.path.join(BB, sec)
    if not os.path.isdir(dst):
        with tarfile.open(os.path.join(D, '..', f'section-{sec}', 'runs', f'runs{sec}.tar.xz')) as tf:
            tf.extractall(dst, members=[m for m in tf.getmembers() if any(s in m.name for s in subs) and ('_ops_' in m.name or '_pps_' in m.name)])
SYS = {'ops': ('F1_ops', 'o-Ps (PHEN=2)'), 'pps': ('F0_pps', 'p-Ps (PHEN=1)')}
txt = lambda p: open(p).read()
same = open(f'{R}/K/D_F1_ops_u5c_s2048.err', 'rb').read() == open(f'{BB}/846/runs/OZ/D_F1_ops_u5c_s2048.err', 'rb').read()
files = lambda m: [f'{R}/{m}/D_{nm}_{k}_s2048.err' for nm, _ in SYS.values() for k in NODES]
lines = all('PAIR START835' in txt(f) for m in ('OZT', 'ZDCT') for f in files(m)) and all('PAIR ORETZDC846 on' in txt(f) for f in files('OZT')) \
    and all('PAIR SODIRAC841 on' in txt(f) and 'PAIR SOCANON843 on' in txt(f) for f in files('ZDCT'))
pr(f"H862b (K = 846 OZ bajt w bajt; START835 w OZT / ZDCT, ORETZDC846 w OZT, SODIRAC841 + SOCANON843 w ZDCT): {'TAK' if same and lines else 'NIE'}  (K {same}, wiersze {lines})")
PATHS = {'OZT': f'{R}/OZT', 'ZDCT': f'{R}/ZDCT', 'OZ': f'{BB}/846/runs/OZ', 'ZDC': f'{BB}/846/runs/ZDC', 'OS': f'{BB}/836/runs/OS', 'ZS': f'{BB}/835/runs/ZS'}
out = {}
for sk, (nm, lab) in SYS.items():
    V = {(m, hf): {k: read(f'{p}/D_{nm}_{k}_s2048.err', half=hf) for k in NODES} for m, p in PATHS.items() for hf in (False, True)}
    nu = {key: nu_of({k: v['om0'] for k, v in d.items()}) for key, d in V.items()}
    f = lambda a, b, hf=False: nu[(a, hf)] - nu[(b, hf)]
    RT, RZ, RO = f('OZT', 'ZDCT'), f('OZ', 'ZDC'), f('OS', 'ZS')
    out[sk] = (RT, RZ, RO, f('ZDCT', 'ZS'))
    pr(f'  {lab}: R_ZDC,T = OZT - ZDCT {units(RT)}')
    pr(f'    R_ZDC 846 = OZ - ZDC {units(RZ)};  R_old = OS - ZS {units(RO)}')
    pr(f'    R_ZDC,T - R_old {(RT - RO)/1e6:+.4f} MHz;  R_ZDC,T - R_ZDC {(RT - RZ)/1e6:+.4f} MHz;  OZT - OS {f("OZT", "OS")/1e6:+.4f}, ZDCT - ZS {f("ZDCT", "ZS")/1e6:+.4f},'
       f' ZDCT - ZDC {f("ZDCT", "ZDC")/1e6:+.4f}, OZT - OZ {f("OZT", "OZ")/1e6:+.4f} MHz')
    pr(f'    polowa okna: R_ZDC,T {f("OZT", "ZDCT", True)/1e6:+.4f}, R_ZDC {f("OZ", "ZDC", True)/1e6:+.4f}, R_old {f("OS", "ZS", True)/1e6:+.4f} MHz')
    pr('    wezly R_ZDC,T / R_old: ' + ', '.join(f"{k} {(V[('OZT', False)][k]['om0'] - V[('ZDCT', False)][k]['om0'])/(V[('OS', False)][k]['om0'] - V[('ZS', False)][k]['om0']):.3f}" for k in NODES)
       + ';  R_ZDC / R_old: ' + ', '.join(f"{k} {(V[('OZ', False)][k]['om0'] - V[('ZDC', False)][k]['om0'])/(V[('OS', False)][k]['om0'] - V[('ZS', False)][k]['om0']):.3f}" for k in NODES))
    pr(f"    apocentra min OZT {min(v['n'] for v in V[('OZT', False)].values())}, ZDCT {min(v['n'] for v in V[('ZDCT', False)].values())}; zle apsydy OZT {sum(v['nbad'] for v in V[('OZT', False)].values())}, ZDCT {sum(v['nbad'] for v in V[('ZDCT', False)].values())}")
    pr('    nachylenie ln T  OZT: ' + ', '.join(f"{V[('OZT', False)][k]['slope']:+.2e}" for k in NODES) + ';  OZ: ' + ', '.join(f"{V[('OZ', False)][k]['slope']:+.2e}" for k in NODES))
RT, RZ, RO, EZ = out['ops']; d = abs(RT - RO)
pr(f"H862s (GLOWNA, start; o-Ps |R_ZDC,T - R_old| = {d/1e6:.4f} MHz; <= 0,015 TAK, >= 0,030 NIE): {'TAK' if d <= 0.015e6 else ('NIE' if d >= 0.030e6 else 'NIEROZSTRZYGNIETE')}")
q = (RT - RO)/EZ
pr(f"H862a (o-Ps |R_ZDC,T - R_ZDC| = {abs(RT - RZ)/1e6:.4f} <= 0,015 i q = {q:+.4f} w [-0,024; -0,008]): {'TAK' if abs(RT - RZ) <= 0.015e6 and -0.024 <= q <= -0.008 else 'NIE'}")
RT, RZ, RO, EZ = out['pps']
pr(f"H862p (p-Ps |R_ZDC,T - R_old| = {abs(RT - RO)/1e6:.4f}, |R_ZDC,T - R_ZDC| = {abs(RT - RZ)/1e6:.4f}; oba <= 0,015): {'TAK' if abs(RT - RO) <= 0.015e6 and abs(RT - RZ) <= 0.015e6 else 'NIE'}")

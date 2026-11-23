# Audit 837: nature of the 835 residual (mode ZT).  Usage: python3 ana837.py -> out/ana837.txt.  Written before the runs;
# gauss, read, nu_of copied from ana835.py.  Z from 830 (u5a,b,d,e) / 828 (u5c); ZT from 834; ZS from 835 archives.
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
SEL = {'830': lambda n: '/D_' in n and n.endswith('_s2048.err') and 'u5c' not in n,
       '828': lambda n: '/D_' in n and n.endswith('_s2048.err') and 'u5c' in n,
       '834': lambda n: '/ZT/' in n, '835': lambda n: '/ZS/' in n}
for sec, f in SEL.items():
    if not os.path.isdir(os.path.join(BS, sec)):
        with tarfile.open(os.path.join(D, '..', f'section-{sec}', 'runs', f'runs{sec}.tar.xz')) as tf:
            tf.extractall(os.path.join(BS, sec), members=[m for m in tf.getmembers() if f(m.name)])
def path(m, nm, k):
    if m == 'Z': return os.path.join(BS, '828' if k == 'u5c' else '830', 'runs', 'Z', f'D_{nm}_{k}_s2048.err')
    if m == 'ZT': return os.path.join(BS, '834', 'runs', 'ZT', f'D_{nm}_{k}_s2048.err')
    if m == 'ZS': return os.path.join(BS, '835', 'runs', 'ZS', f'D_{nm}_{k}_s2048.err')
    return os.path.join(R, m, f'D_{nm}_{k}_s2048.err')
SYS = {'ops': ('F1_ops', 'o-Ps (PHEN=2)'), 'H0': ('F0_H', 'wodor F = 0 (PHEN=1)'), 'H1': ('F1_H', 'wodor F = 1 (PHEN=2)')}
MODES = ('Z', 'ZT', 'ZS', 'X', 'V', 'E1', 'E2', 'LS')
same = open(f'{R}/K/D_F1_ops_u5c_s2048.err', 'rb').read() == open(path('ZS', 'F1_ops', 'u5c'), 'rb').read()
pr(f"H837b (K = ZS 835 bajt w bajt): {'TAK' if same else 'NIE'}")
for nm in ('F1_ops', 'F0_H', 'F1_H'):
    for ln in open(path('E1', nm, 'u5c')):
        if 'START835' in ln: pr(f'  E1 {nm} (tylko czastka 1): {ln.strip()}')
    for ln in open(path('E2', nm, 'u5c')):
        if 'START835' in ln: pr(f'  E2 {nm} (tylko czastka 2): {ln.strip()}')
f = {}; V = {}; okF = True
for sk, (nm, lab) in SYS.items():
    for m in MODES:
        V[(sk, m)] = {k: read(path(m, nm, k)) for k in NODES}
        okF &= all(v['has_sum'] and v['n'] >= 195 and v['nbad'] == 0 for v in V[(sk, m)].values())
    nu = {m: nu_of({k: V[(sk, m)][k]['om0'] for k in NODES}) for m in MODES}
    for m in MODES[1:]: f[(sk, m)] = nu[m] - nu['Z']
    pr(f'  {lab}: ' + '; '.join(f'f({m}) {f[(sk, m)]/1e6:+.4f}' for m in MODES[1:]) + ' MHz')
pr(f'  pliki: wszystkie z linia sum, >= 195 apocentrami i 0 zlych apsyd: {okF}')
okN = True
for sk, (nm, lab) in SYS.items():
    z2 = V[(sk, 'Z')]['u5a']['om0']; s2 = V[(sk, 'ZS')]['u5a']['om0']; t2 = V[(sk, 'ZT')]['u5a']['om0']
    z4 = read(f'{R}/N/Z_{nm}_u5a_s4096.err')['om0']; s4 = read(f'{R}/N/ZS_{nm}_u5a_s4096.err')['om0']
    r2 = s2/z2 - 1; r4 = s4/z4 - 1; okN &= abs(r4 - r2) <= 0.25*abs(r2)
    pr(f'  N {lab} u5a: r(ZS - Z) krok 2048 {r2:+.4e}, 4096 {r4:+.4e} (ZT - Z 2048 {t2/z2 - 1:+.4e})')
pr(f"H837n (|r4096 - r2048| <= 0,25 |r2048|): {'TAK' if okN else 'NIE'}")
okA = True
for sk, (nm, lab) in SYS.items():
    dZ, dS = f[(sk, 'ZT')], f[(sk, 'ZS')]
    ax = f[(sk, 'X')] + f[(sk, 'V')] - dZ; ae = f[(sk, 'E1')] + f[(sk, 'E2')] - dZ
    okA &= abs(ax - dS) <= 0.05*abs(dZ) and abs(ae - dS) <= 0.05*abs(dZ)
    cx, cv = f[(sk, 'X')] - dZ, f[(sk, 'V')] - dZ
    pr(f'  A {lab}: X + V - delta_Z {ax/1e6:+.4f}, E1 + E2 - delta_Z {ae/1e6:+.4f} wobec delta_S {dS/1e6:+.4f} MHz')
    pr(f'    c_x {units(cx)}; c_v {units(cv)}; c_x / c_v {cx/cv:+.4f}; c_x/(c_x+c_v) {cx/(cx+cv):+.4f}')
    pr(f'    lambda* gdyby zla tylko x: {(-dZ - cv)/cx:.4f}; gdyby zla tylko v: {(-dZ - cx)/cv:.4f}')
    pr(f'    c_1 (czastka 1) {units(f[(sk, "E1")] - dZ)}; c_2 (czastka 2) {units(f[(sk, "E2")] - dZ)}')
pr(f"H837a (addytywnosc X + V i E1 + E2 do 0,05 |delta_Z|): {'TAK' if okA else 'NIE'}")
dZ, dS = f[('H0', 'ZT')], f[('H0', 'ZS')]; okW = True
for sk in ('H0', 'H1'):
    dZ, dS = f[(sk, 'ZT')], f[(sk, 'ZS')]
    okW &= abs(f[(sk, 'E2')] - dS) <= 0.25*abs(dS) and abs(f[(sk, 'E1')] - dZ) <= 0.02*abs(dZ)
    pr(f'  W {SYS[sk][1]}: f(E2) {units(f[(sk, "E2")])} wobec delta_S {dS/1e6:+.4f}; f(E1) - delta_Z {units(f[(sk, "E1")] - dZ)}')
pr(f"H837w (GLOWNA: w wodorze reszte niesie czastka 2, przesuniecie czastki 1 nieistotne): {'TAK' if okW else 'NIE'}")
okL = True
for sk, (nm, lab) in SYS.items():
    dS = f[(sk, 'ZS')]; okL &= abs(f[(sk, 'LS')]) <= 0.25*abs(dS)
    pr(f'  L {lab}: f(LS) {units(f[(sk, "LS")])} wobec delta_S {dS/1e6:+.4f}; wezly r(LS - Z) / r(ZT - Z): ' +
       ', '.join(f"{k} {(V[(sk, 'LS')][k]['om0']/V[(sk, 'Z')][k]['om0'] - 1)/(V[(sk, 'ZT')][k]['om0']/V[(sk, 'Z')][k]['om0'] - 1):+.4f}" for k in NODES))
pr(f"H837l (|f(LS)| <= 0,25 |delta_S|): {'TAK' if okL else 'NIE'}")
for m in MODES[1:]:
    P = (f[('H0', m)] + f[('H1', m)])/2; Q = (f[('H1', m)] - f[('H0', m)])/2
    pr(f'  wodor {m}: czesc parzysta {P/1e6:+.4f}, nieparzysta (F1 - F0)/2 {Q/1e6:+.4f} MHz')

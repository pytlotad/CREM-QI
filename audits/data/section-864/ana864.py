# Audit 864: second-order Tulczyjew start (ORDER=2) vs the first-order default (836 OS / 835 ZS archives).  Written before the results.
# Header (units, gauss, read, nu_of) copied from ana863.py.  Usage: LC_ALL=C python3 ana864.py -> out/ana864.txt.
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
import filecmp, shutil
BB = os.path.join(D, 'out', 'base'); os.makedirs(BB, exist_ok=True)
for sec in ('835', '836'):
    with tarfile.open(os.path.join(D, '..', f'section-{sec}', 'runs', f'runs{sec}.tar.xz')) as t: t.extractall(os.path.join(BB, sec), filter='data')
ok = True
for k, ref in (('K_O', '836/runs/OS/D_F1_ops_u5e_s2048'), ('K_Z', '835/runs/ZS/D_F0_H_u5e_s2048'), ('K_1', '836/runs/OS/D_F1_ops_u5e_s2048')):
    s = filecmp.cmp(f'{R}/K/{k}.err', f'{BB}/{ref}.err', shallow=False); ok &= s; pr(f'  {k} = {ref}: {"TAK" if s else "NIE"}')
pr(f'H864b (K_O, K_Z, K_1 bajt w bajt): {"TAK" if ok else "NIE"}')
# H864c from START864 lines
mx_x = 0.0; mx_v = 0.0; nline = 0; terms = []
for M in ('O', 'Z'):
    for f in sorted(os.listdir(f'{R}/{M}')):
        for ln in open(f'{R}/{M}/{f}'):
            if not ln.startswith('PAIR START864'): continue
            w = ln.split(); g = lambda key: float(w[w.index(key) + 1]); nline += 1
            exact = g('(gamma-1)/(gamma+1)'); mx_x = max(mx_x, abs(g('relx') - exact)/max(abs(exact), 1e-300)); mx_v = max(mx_v, abs(g('relv')))
            terms.append((f'{M}/{f[:-4]}', w[w.index('who') + 1], g('gamma-1'), g('relx'), g('relv'), g('fdotTerm'), g('precTerm')))
pr(f'  wiersze START864: {nline}; max |relx - (gamma-1)/(gamma+1)| / (..) = {mx_x:.3e}; max |relv| = {mx_v:.3e}')
for t in terms: pr(f'    {t[0]:28s} who {t[1]} gamma-1 {t[2]:.3e} relx {t[3]:.3e} relv {t[4]:+.3e} fdot/a {t[5]:.2e} prec/a {t[6]:.2e}')
pr(f'H864c: {"TAK" if nline > 0 and mx_x <= 1e-6 and mx_v < 1e-3 else "NIE"}')
CH = (('o-Ps (PHEN=2)', 'D_F1_ops', 0), ('p-Ps (PHEN=1)', 'D_F0_pps', 1), ('wodor F = 0 (PHEN=1)', 'D_F0_H', 0), ('wodor F = 1 (PHEN=2)', 'D_F1_H', 0))
okS = okR = okP = True
for lab, stem, isp in CH:
    om = {}
    for k, p in (('O2', f'{R}/O'), ('Z2', f'{R}/Z'), ('O1', f'{BB}/836/runs/OS'), ('Z1', f'{BB}/835/runs/ZS')):
        om[k] = {n: read(f'{p}/{stem}_{n}_s2048.err')['om0'] for n in NODES}
    nu = {k: nu_of(v) for k, v in om.items()}
    dO, dZ = nu['O2'] - nu['O1'], nu['Z2'] - nu['Z1']; dR = dO - dZ
    pr(f'{lab}:  nu O1 {nu["O1"]/1e6:.4f}, Z1 {nu["Z1"]/1e6:.4f} MHz')
    pr(f'  Delta O (rzad 2 - rzad 1) {units(dO)}'); pr(f'  Delta Z                   {units(dZ)}')
    pr(f'  Delta R                   {units(dR)}   (R1 {(nu["O1"]-nu["Z1"])/1e6:+.4f} -> R2 {(nu["O2"]-nu["Z2"])/1e6:+.4f} MHz)')
    for n, (L, w) in NODES.items():
        pr(f'    {n}: Delta Om O {(om["O2"][n]-om["O1"][n])/(2*math.pi):+.4e} Hz, Z {(om["Z2"][n]-om["Z1"][n])/(2*math.pi):+.4e} Hz (bez wagi wezla)')
    if isp: okP &= abs(dO) <= 10 and abs(dZ) <= 10
    okS &= abs(dO) <= 2e3 and abs(dZ) <= 2e3
    okR &= abs(dR) <= 100
pr(f'H864s (GLOWNA; |Delta| <= 0,002 MHz, 4 kanaly): {"TAK" if okS else "NIE"}')
pr(f'H864r (|Delta R| <= 1e-4 MHz, 4 kanaly): {"TAK" if okR else "NIE"}')
pr(f'H864p (p-Ps |Delta| <= 1e-5 MHz): {"TAK" if okP else "NIE"}')
shutil.rmtree(BB)

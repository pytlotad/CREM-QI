# Audit 867: proper dipole + BMT factor (841) in mode O (OD) and in the conservative Tulczyjew sum (ZD) vs 836 OS / 835 ZS.  Written before
# the results; header copied from ana864.py.  Usage: LC_ALL=C python3 ana867.py -> out/ana867.txt.
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
import filecmp, shutil, re
SIG = 0.015e6   # Hz
BB = os.path.join(D, 'out', 'base'); os.makedirs(BB, exist_ok=True)
for sec in ('835', '836'):
    with tarfile.open(os.path.join(D, '..', f'section-{sec}', 'runs', f'runs{sec}.tar.xz')) as t: t.extractall(os.path.join(BB, sec), filter='data')
ok = filecmp.cmp(f'{R}/K/K_O.err', f'{BB}/836/runs/OS/D_F1_ops_u5e_s2048.err', shallow=False) and \
     filecmp.cmp(f'{R}/K/K_Z.err', f'{BB}/835/runs/ZS/D_F0_H_u5e_s2048.err', shallow=False)
pr(f'H867b (K_O = 836 OS, K_Z = 835 ZS bajt w bajt): {"TAK" if ok else "NIE"}')
nlines = {M: sum(1 for f in os.listdir(f'{R}/{M}') for ln in open(f'{R}/{M}/{f}') if ln.startswith(('PAIR ORETDIRAC867', 'PAIR SODIRAC841'))) for M in ('OD', 'ZD')}
pr(f'  wiersze ORETDIRAC867 w OD / SODIRAC841 w ZD (laczne): {nlines}')
ph = open(os.path.join(D, '..', 'section-840', 'out', 'posthoc840.txt')).read()
blk = re.split(r'odniesienie z 1/2', ph)
g = lambda txt, k: [float(x) for x in re.search(rf'\s{k}\s+F0 ([+-][\d.]+) F1 ([+-][\d.]+)', txt).groups()]
ZT_D, ZT0B_D = g(blk[0], 'ZT'), g(blk[1], 'ZT0')
stat = [(ZT0B_D[i] - ZT_D[i])*1e6 for i in (0, 1)]
pr(f'  statycznie (posthoc840) ZT - D {ZT_D}, ZT0B - D {ZT0B_D} MHz -> Delta_stat F0 {stat[0]/1e6:+.4f}, F1 {stat[1]/1e6:+.4f} MHz')
CH = (('o-Ps (PHEN=2)', 'D_F1_ops'), ('p-Ps (PHEN=1)', 'D_F0_pps'), ('wodor F = 0 (PHEN=1)', 'D_F0_H'), ('wodor F = 1 (PHEN=2)', 'D_F1_H'))
okC = okO = okR = okD = True
for lab, stem in CH:
    om = {k: {n: read(f'{p}/{stem}_{n}_s2048.err')['om0'] for n in NODES} for k, p in
          (('OS', f'{BB}/836/runs/OS'), ('ZS', f'{BB}/835/runs/ZS'), ('OD', f'{R}/OD'), ('ZD', f'{R}/ZD'))}
    nu = {k: nu_of(v) for k, v in om.items()}
    dZ, dO = nu['ZD'] - nu['ZS'], nu['OD'] - nu['OS']; R0, R1 = nu['OS'] - nu['ZS'], nu['OD'] - nu['ZD']
    pr(f'{lab}:')
    pr(f'  Delta_Z (ZD - ZS) {units(dZ)}'); pr(f'  Delta_O (OD - OS) {units(dO)}   Delta_O / Delta_Z {dO/dZ:.4f}')
    pr(f'  R: OS - ZS {R0/1e6:+.4f} -> OD - ZD {R1/1e6:+.4f} MHz, zmiana {units(R1 - R0)}')
    for n in NODES: pr(f'    {n}: Delta_Z {(om["ZD"][n]-om["ZS"][n])/(2*math.pi):+.1f} Hz, Delta_O {(om["OD"][n]-om["OS"][n])/(2*math.pi):+.1f} Hz (bez wagi)')
    okR &= abs(R1 - R0) <= 0.015e6
    if stem == 'D_F0_pps': okO &= abs(dO - dZ) <= 0.005e6
    else: okO &= abs(dO - dZ) <= 0.02*abs(dZ)
    if stem in ('D_F0_H', 'D_F1_H'):
        i = 0 if stem == 'D_F0_H' else 1
        okC &= abs(dZ - stat[i]) <= 0.02*abs(stat[i])
        dev = ZT_D[i]*1e6 + dZ
        pr(f'  Delta_Z / statyczne {dZ/stat[i]:.4f};  odchylka od Diraca po zmianie (ZT - D) + Delta_Z = {units(dev)} = {dev/SIG:+.1f} sigma_2P (bylo {ZT_D[i]*1e6/SIG:+.1f})')
        okD &= (0.05e6 <= dev <= 0.11e6) if i == 0 else (-0.13e6 <= dev <= -0.07e6)
pr(f'H867c (Delta_Z = statyczne do 2 %): {"TAK" if okC else "NIE"}')
pr(f'H867o (GLOWNA; Delta_O = Delta_Z do 2 %, p-Ps do 0,005 MHz): {"TAK" if okO else "NIE"}')
pr(f'H867r (|zmiana R| <= 0,015 MHz): {"TAK" if okR else "NIE"}')
pr(f'H867d (odchylka od Diraca w oknach ZT0B): {"TAK" if okD else "NIE"}')
shutil.rmtree(BB)

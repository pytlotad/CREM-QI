# Audit 865: noise floor of mode O in Ps (start scaled by lambda).  Written before the results; header copied from ana864.py.
# Usage: LC_ALL=C python3 ana865.py -> out/ana865.txt.
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
k = filecmp.cmp(f'{R}/K/K_O.err', f'{BB}/836/runs/OS/D_F1_ops_u5a_s2048.err', shallow=False)
pr(f'H865b (K = 836 OS o-Ps u5a bajt w bajt): {"TAK" if k else "NIE"}')
D864 = {'D_F1_ops': -250.2, 'D_F0_pps': 94.2}   # |Delta R 864| inputs [Hz] from out/ana864.txt (lambda in cm -> Hz), signed
okZ = True; okS = False; dec = True
for lab, stem in (('o-Ps (PHEN=2)', 'D_F1_ops'), ('p-Ps (PHEN=1)', 'D_F0_pps')):
    om = {k: {n: read(f'{p}/{stem}_{n}_s2048.err')['om0'] for n in NODES} for k, p in
          (('O1', f'{BB}/836/runs/OS'), ('Z1', f'{BB}/835/runs/ZS'), ('O_p5', f'{R}/O_p5'), ('O_m5', f'{R}/O_m5'), ('O_p6', f'{R}/O_p6'), ('Z_p5', f'{R}/Z_p5'))}
    d = {k: nu_of({n: om[k][n] - om[('Z1' if k[0] == 'Z' else 'O1')][n] for n in NODES}) for k in ('O_p5', 'O_m5', 'O_p6', 'Z_p5')}
    pr(f'{lab}:')
    for kk in d: pr(f'  Delta {kk:5s} {d[kk]:+10.2f} Hz = {units(d[kk])}')
    for n in NODES:
        pr(f'    {n}: ' + ', '.join(f'{kk} {(om[kk][n] - om[("Z1" if kk[0] == "Z" else "O1")][n])/(2*math.pi):+.2f}' for kk in d) + ' Hz (bez wagi)')
    if stem == 'D_F1_ops': okZ &= -62 <= d['Z_p5'] <= -15
    else: okZ &= abs(d['Z_p5']) <= 1
    c1 = abs(d['O_p5'] - d['Z_p5']); s5 = abs(d['O_p5']) + abs(d['O_m5']); c2 = abs(d['O_p5'] + d['O_m5']) / s5 if s5 else 0.0
    c3 = d['O_p6'] / d['O_p5'] if d['O_p5'] else float('nan')
    i, ii, iii = c1 >= 30, c2 > 0.3, not (0.05 <= c3 <= 0.2)
    okS |= i or ii or iii
    F = max(abs(d['O_p5'] - d['Z_p5']), abs(d['O_m5'] + d['Z_p5']), abs(d['O_p6'] - 0.1*d['Z_p5']))
    pr(f'  (i) |O_p5 - Z_p5| = {c1:.2f} Hz -> {i};  (ii) asymetria {c2:.3f} -> {ii};  (iii) O_p6/O_p5 = {c3:.3f} -> {iii}')
    pr(f'  prog szumu F = {F:.2f} Hz = {units(F)};  |Delta R 864| = {abs(D864[stem]):.1f} Hz, / F = {abs(D864[stem])/F if F else float("inf"):.2f}')
    dec &= abs(D864[stem]) <= 2*F
pr(f'H865z (suma Z liniowa): {"TAK" if okZ else "NIE"}')
pr(f'H865s (GLOWNA; prog szumu O ponad odpowiedz liniowa): {"TAK" if okS else "NIE"}')
pr(f'H865d: {"|Delta R 864| <= 2 F w obu kanalach -> 864 ponizej progu szumu, nie na domyslny" if dec else "skutek 864 ponad 2 F -> fizyczny, kandydat na domyslny"}')
shutil.rmtree(BB)

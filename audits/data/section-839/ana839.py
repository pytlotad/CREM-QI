# Audit 839: analytic D (conjugacy residual, probe p839) against D measured in 838 (ZSA - ZA, apocentre start).  Written before the runs.
import os, math, tarfile, glob
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(D, 'out', 'runs'); BS = os.path.join(D, 'out', 'base')
HEV, C = 4.135667696e-15, 299792458.0
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:+.4f} MHz = {f*HEV*1e6:+.5f} ueV' + (f' (|lambda| {C/abs(f)*100:.4g} cm)' if f else '')
def gauss(N):
    x, w = np.polynomial.legendre.leggauss(N); u = 0.75 + 0.25*x[::-1]; return 1/u, 0.25*w[::-1]
NODES = dict(zip([f'u5{c}' for c in 'abcde'], zip(*gauss(5))))
nu_of = lambda Om: sum(w*L**2*Om[k] for k, (L, w) in NODES.items())/(2*math.pi)
def read(path, tag='mechB0'):   # = ana838.read (om0 only)
    rows = []
    for ln in open(path):
        if not ln.startswith(f'PAIR APSE808 {tag} '): continue
        w = ln.split(); tc = float(w[w.index('tc') + 1])
        r = [np.array([float(v) for v in w[w.index(f'r{q}') + 1:w.index(f'r{q}') + 4]]) for q in range(3)]
        rows.append((tc, r))
    r0 = rows[0][1]; Lh = np.cross(r0[0], r0[2]); Lh /= np.linalg.norm(Lh)
    e1 = r0[1]/np.linalg.norm(r0[1]); e1 -= Lh*np.dot(Lh, e1); e1 /= np.linalg.norm(e1); e2 = np.cross(Lh, e1)
    T = []; P = []; prev = None
    for tc, r in rows:
        th = np.unwrap([math.atan2(np.dot(x, e2), np.dot(x, e1)) for x in r]); rho = [np.linalg.norm(x) for x in r]
        x = th - th[1]; Am = np.vstack([np.ones(3), x, x**2]).T; a0, b, cq = np.linalg.solve(Am, rho)
        phi = th[1] - b/(2*cq)
        if prev is not None: phi += 2*math.pi*round((prev - phi)/(2*math.pi))
        prev = phi; T.append(tc); P.append(phi)
    return np.polyfit(np.array(T), np.array(P), 2)[1]
if not os.path.isdir(os.path.join(BS, 'runs')):
    with tarfile.open(os.path.join(D, '..', 'section-838', 'runs', 'runs838.tar.xz')) as tf:
        tf.extractall(BS, members=[m for m in tf.getmembers() if '/ZA/' in m.name or '/ZSA/' in m.name])
SYS = {'ops': ('F1_ops', 'o-Ps (PHEN=2)'), 'H0': ('F0_H', 'wodor F = 0 (PHEN=1)'), 'H1': ('F1_H', 'wodor F = 1 (PHEN=2)')}
SETS = [0, 1, 3, 7, 15]; SETN = {0: 'P0 Coulomb + Thomas', 1: 'P1 + gamma', 3: 'P2 + Darwin', 7: 'P3 + sily spinowe', 15: 'P4 + precesja (pelny)'}
def om(sk, nd, var):
    out = {}; rr = {}
    for ln in open(f'{R}/{sk}_{nd}_{var}.txt'):
        if ln.startswith('OMR'):
            w = ln.split(); s = int(w[2]); out[s] = -float(w[4]); rr[s] = float(w[6])/float(w[8])
    return out, rr
okC = okR = okN = okD = okW = True; Dan = {}; Dme = {}
for f in sorted(glob.glob(f'{R}/*_CL.txt')):
    for ln in open(f):
        if ln.startswith('CHECK'): okC &= float(ln.split()[-1]) <= 1e-12
pr(f"H839c (suma sondy = allExternalForces, forma L i T, <= 1e-12): {'TAK' if okC else 'NIE'}")
for sk, (nm, lab) in SYS.items():
    base = {nd: om(sk, nd, 'base') for nd in NODES}
    dme = {nd: read(f'{BS}/runs/ZSA/D_{nm}_{nd}_s2048.err') - read(f'{BS}/runs/ZA/D_{nm}_{nd}_s2048.err') for nd in NODES}
    Dme[sk] = nu_of(dme)
    Dan[sk] = {s: nu_of({nd: base[nd][0][s] for nd in NODES}) for s in SETS}
    nN = nu_of({nd: om(sk, nd, 'N2')[0][15] for nd in NODES}); nh = nu_of({nd: om(sk, nd, 'h2')[0][15] for nd in NODES})
    rmax = max(base[nd][1][15] for nd in NODES); okR &= rmax <= 0.01
    okN &= abs(nN - Dan[sk][15]) <= 0.01*abs(Dan[sk][15]) and abs(nh - Dan[sk][15]) <= 0.01*abs(Dan[sk][15])
    okD &= abs(Dan[sk][15] - Dme[sk]) <= 0.15*abs(Dme[sk])
    pr(f'  {lab}:')
    pr(f'    D_an (P4) {units(Dan[sk][15])}   D_838 {units(Dme[sk])}   D_an / D_838 {Dan[sk][15]/Dme[sk]:.4f}')
    pr(f'    numeryka: N 8192 {nN/1e6:+.5f}, h/2 {nh/1e6:+.5f} MHz;  max rms|R| / rms|a^T - a^L| {rmax:.2e}')
    for nd in NODES:
        q = base[nd][0][15]/dme[nd]; okW &= abs(q - 1) <= 0.25
        pr(f'    {nd}: D_an {base[nd][0][15]:+.4e} rad/s  D_838 {dme[nd]:+.4e} rad/s  stosunek {q:.4f}   (' +
           ', '.join(f'P{i} {base[nd][0][s]:+.3e}' for i, s in enumerate(SETS)) + ')')
    prev = 0.0
    for s in SETS:
        pr(f'    {SETN[s]:24s} nu {Dan[sk][s]/1e6:+.4f} MHz  przyrost {units(Dan[sk][s] - prev)}'); prev = Dan[sk][s]
pr(f"H839r (rms|R| <= 0,01 rms|a^T - a^L|): {'TAK' if okR else 'NIE'}")
pr(f"H839n (N x 2, h / 2: <= 1 %): {'TAK' if okN else 'NIE'}")
pr(f"H839d (GLOWNA: |D_an - D_838| <= 0,15 |D_838|): {'TAK' if okD else 'NIE'}")
pr(f"H839w (wezly <= 25 %): {'TAK' if okW else 'NIE'}")
Pan = (Dan['H0'][15] + Dan['H1'][15])/2; Pme = (Dme['H0'] + Dme['H1'])/2
pr(f"H839p (czesc parzysta wodoru {Pan/1e6:+.4f} wobec {Pme/1e6:+.4f} (838) i 0,0630; |P_an - 0,0630| <= 0,25 x 0,0630): " +
   ('TAK' if abs(Pan/1e6 - 0.0630) <= 0.25*0.0630 else 'NIE'))
for s in SETS:
    pr(f'  wodor {SETN[s]:24s}: parzysta {(Dan["H0"][s] + Dan["H1"][s])/2/1e6:+.4f}, nieparzysta (F1 - F0)/2 {(Dan["H1"][s] - Dan["H0"][s])/2/1e6:+.4f} MHz')
pr(f'  838 zmierzone: parzysta {Pme/1e6:+.4f}, nieparzysta {(Dme["H1"] - Dme["H0"])/2/1e6:+.4f} MHz')

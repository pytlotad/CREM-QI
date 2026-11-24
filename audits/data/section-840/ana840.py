# Audit 840: engine Thomas forms vs the classical limit of the Dirac FW Hamiltonian (hydrogen test particle).  Written before the runs.
import os, math, glob
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(D, 'out', 'runs')
HEV, C = 4.135667696e-15, 299792458.0
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:+.4f} MHz = {f*HEV*1e6:+.5f} ueV' + (f' (|lambda| {C/abs(f)*100:.4g} cm)' if f else '')
def gauss(N):
    x, w = np.polynomial.legendre.leggauss(N); u = 0.75 + 0.25*x[::-1]; return 1/u, 0.25*w[::-1]
NODES = dict(zip([f'u5{c}' for c in 'abcde'], zip(*gauss(5))))
nu_of = lambda Om: sum(w*L**2*Om[k] for k, (L, w) in NODES.items())/(2*math.pi)
FORMS = ['Z', 'ZT', 'Z0', 'ZT0']
def rd(sk, nd, var):
    o = {}; rr = {}; sc = None
    for ln in open(f'{R}/{sk}_{nd}_{var}.txt'):
        w = ln.split()
        if ln.startswith('SCALE'): sc = (float(w[2]), float(w[4]))
        if ln.startswith('OMX'): o[w[1]] = -float(w[3]); rr[w[1]] = float(w[5])
    return o, rr, sc
okC = all(float(l.split()[-1]) <= 1e-12 for f in glob.glob(f'{R}/*_CL.txt') for l in open(f) if l.startswith('CHECK'))
pr(f"H840c (suma sondy = allExternalForces, L i T, <= 1e-12): {'TAK' if okC else 'NIE'}")
okR = okN = True; nu = {}; NM = {'H0': 'wodor F = 0 (PHEN=1)', 'H1': 'wodor F = 1 (PHEN=2)'}
for sk in ('H0', 'H1'):
    B = {nd: rd(sk, nd, 'base') for nd in NODES}; N2 = {nd: rd(sk, nd, 'N2') for nd in NODES}; H2 = {nd: rd(sk, nd, 'h2') for nd in NODES}
    pr(f'  {NM[sk]} (czastka probna):')
    for nd in NODES:
        okR &= all(B[nd][1][x] <= 0.01*B[nd][2][0] for x in FORMS)
        pr(f'    {nd}: rms first {B[nd][2][0]:.3e}, spinless |a^D - a^Z| {B[nd][2][1]:.3e}; ' + ', '.join(f'{x} {B[nd][0][x]:+.4e} rad/s (rmsR/first {B[nd][1][x]/B[nd][2][0]:.1e})' for x in FORMS))
    for x in FORMS:
        nu[sk, x] = nu_of({nd: B[nd][0][x] for nd in NODES})
        a2 = nu_of({nd: N2[nd][0][x] for nd in NODES}); b2 = nu_of({nd: H2[nd][0][x] for nd in NODES})
        tol = max(0.01*abs(nu[sk, x]), 1e3); okN &= abs(a2 - nu[sk, x]) <= tol and abs(b2 - nu[sk, x]) <= tol
        pr(f'    nu_{x} - nu_D {units(nu[sk, x])}   (N x 2 {a2/1e6:+.5f}, h / 2 {b2/1e6:+.5f})')
pr(f"H840r (rms|R_X| <= 0,01 pierwszego rzedu): {'TAK' if okR else 'NIE'}")
pr(f"H840n (N x 2, h / 2): {'TAK' if okN else 'NIE'}")
D839 = {'H0': 0.3583e6, 'H1': -0.2306e6}
okK = all(abs((nu[sk, 'ZT'] - nu[sk, 'Z']) - D839[sk]) <= 0.03*abs(D839[sk]) for sk in D839)
pr('  ZT - Z (czastka probna) ' + ', '.join(f'{sk} {(nu[sk, "ZT"] - nu[sk, "Z"])/1e6:+.4f} (839 {D839[sk]/1e6:+.4f})' for sk in D839) + ' MHz')
pr(f"H840k (spojnosc z 839, 3 %): {'TAK' if okK else 'NIE'}")
pr('H840m (GLOWNA: |nu_X - nu_D| <= 0,015 MHz w F = 0 i F = 1):')
for x in FORMS:
    ok = all(abs(nu[sk, x]) <= 0.015e6 for sk in ('H0', 'H1'))
    pr(f'  {x:4s}: F = 0 {nu["H0", x]/1e6:+.4f}, F = 1 {nu["H1", x]/1e6:+.4f} MHz; parzysta {(nu["H0", x] + nu["H1", x])/2e6:+.4f}, '
       f'nieparzysta (F1 - F0)/2 {(nu["H1", x] - nu["H0", x])/2e6:+.4f} MHz -> {"TAK" if ok else "NIE"}')

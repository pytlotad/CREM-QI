# Audit 843: CREM_SPIN_ORBIT_CANONICAL on top of 841 -- static S (vs Dirac FW), static P (two-body Delta_C), dynamic p814 ZDCA - ZDA(841).
# Written before the runs (from ana841.py).
import os, math, tarfile
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(D, 'out', 'runs'); BS = os.path.join(D, 'out', 'base')
HEV, C = 4.135667696e-15, 299792458.0
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:+.4f} MHz = {f*HEV*1e6:+.5f} ueV' + (f' (|lambda| {C/abs(f)*100:.4g} cm)' if f else '')
def gauss(N):
    x, w = np.polynomial.legendre.leggauss(N); u = 0.75 + 0.25*x[::-1]; return 1/u, 0.25*w[::-1]
NODES = dict(zip([f'u5{c}' for c in 'abcde'], zip(*gauss(5))))
nu_of = lambda Om: sum(w*L**2*Om[k] for k, (L, w) in NODES.items())/(2*math.pi)
def read(path, tag='mechB0'):   # = ana839.read
    rows = []
    for ln in open(path):
        if not ln.startswith(f'PAIR APSE808 {tag} '): continue
        w = ln.split(); tc = float(w[w.index('tc') + 1])
        r = [np.array([float(v) for v in w[w.index(f'r{q}') + 1:w.index(f'r{q}') + 4]]) for q in range(3)]
        rows.append((tc, r))
    r0 = rows[0][1]; Lh = np.cross(r0[0], r0[2]); Lh /= np.linalg.norm(Lh)
    e1 = r0[1]/np.linalg.norm(r0[1]); e1 -= Lh*np.dot(Lh, e1); e1 /= np.linalg.norm(e1); e2 = np.cross(Lh, e1)
    T = []; P = []; prev = None; nbad = 0
    for tc, r in rows:
        th = np.unwrap([math.atan2(np.dot(x, e2), np.dot(x, e1)) for x in r]); rho = [np.linalg.norm(x) for x in r]
        x = th - th[1]; Am = np.vstack([np.ones(3), x, x**2]).T; a0, b, cq = np.linalg.solve(Am, rho)
        xs = -b/(2*cq); hh = max(abs(x[0]), abs(x[2]))
        if not (cq < 0 and x[0] - 0.5*hh <= xs <= x[2] + 0.5*hh): nbad += 1
        phi = th[1] + xs
        if prev is not None: phi += 2*math.pi*round((prev - phi)/(2*math.pi))
        prev = phi; T.append(tc); P.append(phi)
    return np.polyfit(np.array(T), np.array(P), 2)[1], len(T), nbad
if not os.path.isdir(os.path.join(BS, 'runs')):
    with tarfile.open(os.path.join(D, '..', 'section-841', 'runs', 'runs841.tar.xz')) as tf:
        tf.extractall(BS, members=[m for m in tf.getmembers() if '/ZDA/' in m.name])
same = open(f'{R}/K/D_F0_H_u5c_s2048.err', 'rb').read() == open(f'{BS}/runs/ZDA/D_F0_H_u5c_s2048.err', 'rb').read()
printed = all('PAIR SOCANON843 on' in open(f'{R}/ZDCA/D_{nm}_{nd}_s2048.err').read() for nm in ('F1_ops', 'F0_H', 'F1_H') for nd in NODES)
nosod = 'SOCANON843' not in open(f'{R}/K/D_F0_H_u5c_s2048.err').read()
pr(f"H843b (K = ZDA 841 bajt w bajt; SOCANON843 we wszystkich ZDCA, brak w K): {'TAK' if same and printed and nosod else 'NIE'}")
# static S
def omx(f):
    o = {}
    for ln in open(f):
        if ln.startswith('OMX'): w = ln.split(); o[w[1]] = -float(w[3])
    return o
S = {}
for sk in ('H0', 'H1'):
    for var in ('', '_h2'):
        O = {nd: omx(f'{R}/S/{sk}_{nd}{var}.txt') for nd in NODES}
        for x in ('ZD', 'ZDC'): S[sk, var, x] = nu_of({nd: O[nd][x] for nd in NODES})
for x in ('ZD', 'ZDC'):
    pr(f'  S: nu_{x} - nu_D: F = 0 {units(S["H0", "", x])}, F = 1 {units(S["H1", "", x])};  parzysta {(S["H0", "", x] + S["H1", "", x])/2e6:+.4f}, '
       f'nieparzysta {(S["H1", "", x] - S["H0", "", x])/2e6:+.4f} MHz  (h/2: {S["H0", "_h2", x]/1e6:+.4f} / {S["H1", "_h2", x]/1e6:+.4f})')
okS = all(abs(S[sk, '', 'ZDC']) <= 0.015e6 for sk in ('H0', 'H1'))
pr(f"H843s (GLOWNA: |nu_ZDC - nu_D| <= 0,015 MHz w F = 0 i F = 1): {'TAK' if okS else 'NIE'}  (F = 0 {S['H0', '', 'ZDC']/1e6:+.5f}, F = 1 {S['H1', '', 'ZDC']/1e6:+.5f} MHz)")
# static P and dynamic M
D841 = {'ops': -0.0745e6, 'H0': -0.5261e6, 'H1': +0.5260e6}   # 841: nu(ZDA) - nu(ZA), p814
SYS = {'ops': ('F1_ops', 'o-Ps (PHEN=2)'), 'H0': ('F0_H', 'wodor F = 0 (PHEN=1)'), 'H1': ('F1_H', 'wodor F = 1 (PHEN=2)')}
okM = True
for sk, (nm, lab) in SYS.items():
    P = {}; P2 = {}
    for nd in NODES:
        w = open(f'{R}/P/{sk}_{nd}.txt').read().split(); P[nd] = {q: float(w[w.index(q) + 1]) for q in ('ZD',)}
        w = open(f'{R}/P/{sk}_{nd}_N2.txt').read().split(); P2[nd] = float(w[w.index('ZD') + 1])
    dP = {q: nu_of({nd: P[nd][q] for nd in NODES}) for q in ('ZD',)}; dP2 = nu_of(P2)
    Om = {}; bad = 0; nmin = 999
    for nd in NODES:
        z, n1, b1 = read(f'{R}/ZDCA/D_{nm}_{nd}_s2048.err'); a, n2, b2 = read(f'{BS}/runs/ZDA/D_{nm}_{nd}_s2048.err')
        Om[nd] = z - a; bad += b1 + b2; nmin = min(nmin, n1, n2)
    dM = nu_of(Om); okM &= abs(dM - dP['ZD']) <= 0.15*abs(dP['ZD'])
    pr(f'  {lab}: Delta_C statycznie (ZDC - ZD) {units(dP["ZD"])} [N x 2 {dP2/1e6:+.5f}]; lacznie ZDC - Z (841 + 843) {units(D841[sk] + dM)}')
    pr(f'    Delta_C p814 {units(dM)}  stosunek {dM/dP["ZD"]:.4f}; apocentra min {nmin}, zle apsydy {bad}')
    pr('    wezly p814 / P: ' + ', '.join(f'{nd} {Om[nd]/P[nd]["ZD"]:.3f}' for nd in NODES))
pr(f"H843m (|Delta_C p814 - Delta_C stat.| <= 0,15 |Delta_C stat.|): {'TAK' if okM else 'NIE'}")

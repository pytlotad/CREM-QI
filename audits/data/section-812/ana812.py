# Audit 812: Ps 2P orbit torus difference E_orb(2,2) - E_orb(2,1) from the apsidal precession integrated over L (method 808) ->
# <r^-3> of the Ps triplet as a measured orbit property instead of Born's rule.  Usage: python3 ana812.py [runs dir] (default out/runs)
import os, sys, math, filecmp, tarfile, io
import numpy as np
D = os.path.dirname(os.path.abspath(__file__))
O = sys.argv[1] if len(sys.argv) > 1 else os.path.join(D, 'out', 'runs')
HEV, C, h = 4.135667696e-15, 299792458.0, 6.62607015e-34
NU4 = 5474.642e6                                   # m_e c^2 alpha^4 / 64 (811)
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.4f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
def gauss(N):
    x, w = np.polynomial.legendre.leggauss(N); u = 0.75 + 0.25*x[::-1]; return 1/u, 0.25*w[::-1]
N3 = dict(zip(('u3a', 'u3b', 'u3c'), zip(*gauss(3)))); N4 = dict(zip(('u4a', 'u4b', 'u4c', 'u4d'), zip(*gauss(4))))
def apse(path, tag='mechB0'):   # = ana808.apse
    rows = []; has_sum = False
    for ln in open(path):
        if ln.startswith(f'PAIR TORQUE799 {tag} sum'): has_sum = True
        if not ln.startswith(f'PAIR APSE808 {tag} '): continue
        w = ln.split(); tc = float(w[w.index('tc') + 1])
        r = [np.array([float(v) for v in w[w.index(f'r{q}') + 1:w.index(f'r{q}') + 4]]) for q in range(3)]
        rows.append((tc, r))
    if not rows: return float('nan'), 0, 0, float('nan'), has_sum
    r0 = rows[0][1]; Lh = np.cross(r0[0], r0[2]); Lh /= np.linalg.norm(Lh)
    e1 = r0[1]/np.linalg.norm(r0[1]); e1 -= Lh*np.dot(Lh, e1); e1 /= np.linalg.norm(e1); e2 = np.cross(Lh, e1)
    T = []; P = []; nbad = 0; prev = None
    for tc, r in rows:
        th = np.unwrap([math.atan2(np.dot(x, e2), np.dot(x, e1)) for x in r]); rho = [np.linalg.norm(x) for x in r]
        x = th - th[1]; A = np.vstack([np.ones(3), x, x**2]).T; a0, b, cq = np.linalg.solve(A, rho)
        xs = -b/(2*cq); hh = max(abs(x[0]), abs(x[2]))
        if not (cq < 0 and x[0] - 0.5*hh <= xs <= x[2] + 0.5*hh): nbad += 1
        phi = th[1] + xs
        if prev is not None: phi += 2*math.pi*round((prev - phi)/(2*math.pi))
        prev = phi; T.append(tc); P.append(phi)
    T = np.array(T); P = np.array(P); k, b0 = np.polyfit(T, P, 1); res = P - (k*T + b0)
    return k, len(T), nbad, float(np.sqrt(np.mean(res**2))), has_sum
# --- H812k
ref = os.path.join(D, '..', 'section-808', 'runs', 'runs808.tar.xz'); okK = False
with tarfile.open(ref) as t:
    m = [x for x in t.getnames() if x.endswith('O1/F1_proton_u3b_s2048_5.err')]
    if m: okK = t.extractfile(m[0]).read() == open(f'{O}/K/F1_proton_u3b_s2048_5.err', 'rb').read()
pr(f"H812k (bez przelacznika plik wodoru 808 O1/F1_proton_u3b_s2048_5 bajt w bajt): {'TAK' if okK else 'NIE'}")
# --- H812v
okV = True; nf = 0; minapo = 10**9; nbadv = 0
def AP(path):
    global okV, nf, minapo, nbadv
    k, n, nb, rms, hs = apse(path); nf += 1; minapo = min(minapo, n); nbadv += nb; okV &= hs and n >= 195 and nb == 0; return k
H = (1024, 2048, 4096)
O3 = {(s, k): AP(f'{O}/P3/F1_ops_{k}_s{s}_5.err') for s in H for k in N3}
O4 = {(s, k): AP(f'{O}/P4/F1_ops_{k}_s{s}_5.err') for s in H for k in N4}
PP = {k: AP(f'{O}/PP/F0_pps_{k}_s2048_5.err') for k in N3}
okV &= nf == 21
pr(f"H812v (21 plikow Ps: linia mechB0 sum, >= 195 apocentrow, 0 zlych wierzcholkow): {'TAK' if okV else 'NIE'} (plikow {nf}, apocentrow min {minapo}, zlych {nbadv})")
# --- H812c
dc = max(abs(PP[k]/O3[(2048, k)] - 1) for k in N3)
pr(f"H812c (p-Ps = o-Ps przy 2048 do 1e-9 w kazdym wezle): {'TAK' if dc <= 1e-9 else 'NIE'} (max {dc:.2e})")
for k, (L, w) in N3.items(): pr(f'  {k} L = {L:.6f}: o-Ps {O3[(2048, k)]:.10e}, p-Ps {PP[k]:.10e} rad/s')
# --- H812s
nu_of = lambda nodes, Om: sum(w*L**2*Om[k] for k, (L, w) in nodes.items())/(2*math.pi)
pobs = lambda y: math.log2((y[0] - y[1])/(y[1] - y[2])) if (y[0] - y[1])*(y[1] - y[2]) > 0 else float('nan')
okS = True; R = {}
for nodes, Om, nm in ((N3, O3, '3'), (N4, O4, '4')):
    R2 = {}; R4 = {}
    for k in nodes:
        y = [Om[(s, k)] for s in H]; p = pobs(y); okS &= 1.9 <= p <= 2.1
        R2[k] = y[2] + (y[2] - y[1])/3; r2p = y[1] + (y[1] - y[0])/3; R4[k] = R2[k] + (R2[k] - r2p)/15
        pr(f'  {nm} wezly {k}: Omega 1024 / 2048 / 4096 = ' + ' / '.join(f'{v:.10e}' for v in y) + f'; p = {p:.4f}')
    ys = [nu_of(nodes, {k: Om[(s, k)] for k in nodes}) for s in H]; n2, n4 = nu_of(nodes, R2), nu_of(nodes, R4)
    okS &= abs(n4 - n2) <= 0.01e6; R[nm] = (ys, n2, n4)
    pr(f'  nu_orb {nm} wezly 1024 / 2048 / 4096: ' + ' / '.join(f'{v/1e6:.4f}' for v in ys) + f' MHz; R2 {n2/1e6:.4f}, R4 {n4/1e6:.4f} MHz')
pr(f"H812s (p w [1,9; 2,1] w kazdym wezle, |R4 - R2| <= 0,01 MHz): {'TAK' if okS else 'NIE'}")
dq = R['3'][2] - R['4'][2]
pr(f"H812q (|3 wezly - 4 wezly| w granicy kroku <= 0,010 MHz): {'TAK' if abs(dq) <= 0.010e6 else 'NIE'} ({dq/1e6:+.4f} MHz)")
# --- H812o
nuO = R['4'][2]; c_o = R['3'][2] - R['3'][0][1]; nuP = nu_of(N3, PP) + c_o
rhoO, rhoP = nuO/NU4, nuP/NU4
for nm, nu, rho in (('o-Ps (PHEN=2, granica kroku, 4 wezly)', nuO, rhoO), ('p-Ps (PHEN=1, 2048 + poprawka kroku o-Ps)', nuP, rhoP)):
    pr(f'  {nm}: nu_orb = {units(nu)}; / (m_e c^2 alfa^4/64) = {rho:.7f}; delta_6 = {(nu - NU4)/1e6:+.4f} MHz (jednocialowy Sommerfeld +0,182)')
okO = abs(rhoO - 1) <= 2e-4 and abs(rhoP - 1) <= 2e-4
pr(f"H812o (GLOWNA; |nu_orb/5 474,642 - 1| <= 2e-4, o-Ps i p-Ps): {'TAK' if okO else 'NIE'}; rho = <r^-3>_orbita/<r^-3>_Born = {rhoO:.7f} "
   f"(klasyczne 8/9 = 0,8888889 -- odleglosc {abs(rhoO - 8/9)/2e-4:.0f} progow)")
# --- H812m: o-Ps multiplet with <r^-3> from the orbit
from mpmath import mp, mpf
mp.dps = 30
me = mpf('9.1093837139e-31'); al = 1/mpf('137.035999177'); cc = mpf(C); hb = mpf(h)/(2*mp.pi); k_ = al*hb*cc
g = mpf('2.00231930436092'); muB = mpf('9.2740100783e-24'); mom = g*muB/2
aPs = hb**2/((me/2)*k_); r3 = 1/(24*aPs**3)
aSO = (2*g - 1)*k_*hb**2*r3/(2*me**2*cc**2); bT = mpf('1e-7')*mom*mom*r3
E = {}
for J in (0, 1, 2):
    LS = mpf(J*(J + 1) - 2 - 2)/2; Ll = mpf(3)/2; L2 = 2
    P = mpf(2)/3 - (LS*(LS + mpf(1)/2) - mpf(2)*L2/3)/(2*(Ll - 1)*(Ll + 1)); E[J] = float((aSO*LS + bT*(6*P - 4))/mpf(h))
ENG = {0: 18496.12e6, 1: 13010.90e6, 2: 8625.26e6}
MEAS = {0: (18501.02e6, 0.61e6, 'Gurung i in., PRL 125, 073002 (2020)'), 1: (13012.42e6, 1.67e6, 'Hagena i in., PRL 71, 2887 (1993)'),
        2: (8627.94e6, 0.95e6, 'Sheldon i in., PRL 131, 043001 (2023)')}
okM = True
pr('o-Ps 2^3S_1 - 2^3P_J z <r^-3> z orbity (rho o-Ps):')
for J in (0, 1, 2):
    d = -(rhoO - 1)*E[J]; okM &= abs(d) <= 0.5e6; nu = ENG[J] + d; m, s, ref_ = MEAS[J]
    pr(f'  J = {J}: E_J (czesc spinowa) {E[J]/1e6:+.3f} MHz; delta nu = {d/1e6:+.4f} MHz; nu_{J} = {units(nu)}; pomiar {m/1e6:.2f}({s/1e6:.2f}) MHz ({ref_}); '
       f'reszta {(nu - m)/1e6:+.2f} MHz = {(nu - m)*HEV*1e6:+.4f} ueV ({(nu - m)/s:+.1f} sigma); przed {(ENG[J] - m)/1e6:+.2f} MHz')
pr(f"H812m (|delta nu_J| <= 0,5 MHz, J = 0, 1, 2): {'TAK' if okM else 'NIE'}")
pr('p-Ps 2^1P_1: bez czesci spinowej -- nic do zmiany (orbita jak o-Ps: H812c).')
allok = okK and okV and dc <= 1e-9 and okS and abs(dq) <= 0.010e6 and okO
pr(f"DECYZJA: {'<r^-3> trypletu Ps = zmierzona roznica torusow orbity, nie regula Borna (przeklasyfikowanie w Model.md)' if allok else 'regula Borna zostaje importem'}")

# Audit 815: the Plummer separation floor as the source of the 814 engine - integral gap, and the eps -> 0 limit of the Ps torus
# difference under the production force.  Written during the runs, before reading them.  Usage: python3 ana815.py [runs dir] -> out/ana815.txt
import os, sys, math, tarfile
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D); sys.path.insert(0, os.path.join(D, '..', 'section-814'))
import pred815 as P15
from pred814 import parse_start
O = sys.argv[1] if len(sys.argv) > 1 else os.path.join(D, 'out', 'runs')
HEV, C = 4.135667696e-15, 299792458.0; al = 1/137.035999177; NU4 = 5474.642e6
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.4f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
x, w = np.polynomial.legendre.leggauss(5); u = 0.75 + 0.25*x[::-1]
NODES = dict(zip([f'u5{c}' for c in 'abcde'], zip(1/u, 0.25*w[::-1])))
nu_of = lambda Om: sum(w*L**2*Om[k] for k, (L, w) in NODES.items())/(2*math.pi)
def read(path, tag='mechB0'):   # = ana814.read (apse vertices as ana808)
    rows = []; has_sum = False; has_cons = False
    for ln in open(path):
        if ln.startswith(f'PAIR TORQUE799 {tag} sum'): has_sum = True
        if ln.startswith(f'PAIR CONS814 {tag} start'): has_cons = True
        if not ln.startswith(f'PAIR APSE808 {tag} '): continue
        w_ = ln.split(); tc = float(w_[w_.index('tc') + 1])
        r = [np.array([float(v) for v in w_[w_.index(f'r{q}') + 1:w_.index(f'r{q}') + 4]]) for q in range(3)]
        rows.append((tc, r))
    r0 = rows[0][1]; Lh = np.cross(r0[0], r0[2]); Lh /= np.linalg.norm(Lh)
    e1 = r0[1]/np.linalg.norm(r0[1]); e1 -= Lh*np.dot(Lh, e1); e1 /= np.linalg.norm(e1); e2 = np.cross(Lh, e1)
    T = []; Pp = []; nbad = 0; prev = None
    for tc, r in rows:
        th = np.unwrap([math.atan2(np.dot(v, e2), np.dot(v, e1)) for v in r]); rho = [np.linalg.norm(v) for v in r]
        xx = th - th[1]; Am = np.vstack([np.ones(3), xx, xx**2]).T; a0, b, cq = np.linalg.solve(Am, rho)
        xs = -b/(2*cq); hh = max(abs(xx[0]), abs(xx[2]))
        if not (cq < 0 and xx[0] - 0.5*hh <= xs <= xx[2] + 0.5*hh): nbad += 1
        phi = th[1] + xs
        if prev is not None: phi += 2*math.pi*round((prev - phi)/(2*math.pi))
        prev = phi; T.append(tc); Pp.append(phi)
    T = np.array(T); Pp = np.array(Pp)
    return dict(lin=np.polyfit(T, Pp, 1)[0], n=len(T), nbad=nbad, has_sum=has_sum, has_cons=has_cons, path=path)
def tarbytes(tarpath, suffix):
    with tarfile.open(tarpath) as t:
        m = [q for q in t.getnames() if q.endswith(suffix)]; return t.extractfile(m[0]).read() if m else None
T814 = os.path.join(D, '..', 'section-814', 'runs', 'runs814.tar.xz')
# ---------------- H815k
k1 = tarbytes(T814, 'P5/F1_ops_u5a_s2048_5.err') == open(f'{O}/K/F1_ops_u5a_s2048_5.err', 'rb').read()
pr(f"H815k (skala 0,05 jawnie = 814 u5a 2048 bajt w bajt): {'TAK' if k1 else 'NIE'}")
# ---------------- load / H815v
A = {(0.025, s, k): read(f'{O}/S025/F1_ops_{k}_s{s}_5.err') for s in (2048, 4096, 8192) for k in NODES}
A.update({(0.1, s, k): read(f'{O}/S100/F1_ops_{k}_s{s}_5.err') for s in (2048, 4096) for k in NODES})
PP = {k: read(f'{O}/PP025/F0_pps_{k}_s2048_5.err') for k in NODES}
allf = list(A.values()) + list(PP.values())
okV = len(allf) == 30 and all(v['has_sum'] and v['has_cons'] and v['n'] >= 195 and v['nbad'] == 0 for v in allf)
pr(f"H815v (30 plikow: sum, CONS814, >= 195 apocentrow, 0 zlych): {'TAK' if okV else 'NIE'} (plikow {len(allf)}, "
   f"min apocentrow {min(v['n'] for v in allf)}, zlych {sum(v['nbad'] for v in allf)})")
# ---------------- H815s (0,025: 2048 / 4096 / 8192)
pobs = lambda y: math.log2((y[0] - y[1])/(y[1] - y[2])) if (y[0] - y[1])*(y[1] - y[2]) > 0 else float('nan')
okS = True; R2a, R4a, R2b = {}, {}, {}
for k in NODES:
    y = [A[(0.025, s, k)]['lin'] for s in (2048, 4096, 8192)]; p = pobs(y); okS &= 1.8 <= p <= 2.2
    R2a[k] = y[2] + (y[2] - y[1])/3; r2p = y[1] + (y[1] - y[0])/3; R4a[k] = R2a[k] + (R2a[k] - r2p)/15
    z = [A[(0.1, s, k)]['lin'] for s in (2048, 4096)]; R2b[k] = z[1] + (z[1] - z[0])/3
    pr(f'  {k}: skala 0,025 Omega 2048/4096/8192 = ' + ' / '.join(f'{v:.10e}' for v in y) + f'; p = {p:.4f}; R2 {R2a[k]:.10e}, R4 {R4a[k]:.10e}')
    pr(f'       skala 0,1 Omega 2048/4096 = ' + ' / '.join(f'{v:.10e}' for v in z) + f'; R2 {R2b[k]:.10e}')
n2a, n4a, n2b = nu_of(R2a), nu_of(R4a), nu_of(R2b); okS &= abs(n4a - n2a) <= 0.010e6
pr(f'  nu_orb skala 0,025: R2 {n2a/1e6:.4f}, R4 {n4a/1e6:.4f} MHz; skala 0,1: R2 {n2b/1e6:.4f} MHz')
pr(f"H815s (0,025: p w [1,8; 2,2], |R4 - R2| <= 0,010 MHz): {'TAK' if okS else 'NIE'} (R4 - R2 {(n4a - n2a)/1e6:+.4f} MHz)")
# ---------------- predictions
okP = True; conv = 0.0; Pa, Pb, P0 = {}, {}, {}
for k in NODES:
    sa = parse_start(A[(0.025, 8192, k)]['path']); sb = parse_start(A[(0.1, 4096, k)]['path'])
    dst = max(np.max(np.abs(sa[q] - sb[q])) for q in ('x1', 'x2', 'v1', 'v2'))
    a12 = P15.omega(sa, 0.025, rtol=1e-12)[0]; Pa[k], nap = P15.omega(sa, 0.025); conv = max(conv, abs(Pa[k]/a12 - 1))
    b12 = P15.omega(sb, 0.1, rtol=1e-12)[0]; Pb[k] = P15.omega(sb, 0.1)[0]; conv = max(conv, abs(Pb[k]/b12 - 1))
    P0[k] = P15.omega(sa, 0.0)[0]
    ra, rb = R4a[k]/Pa[k] - 1, R2b[k]/Pb[k] - 1; okP &= abs(ra) <= 5e-6 and abs(rb) <= 5e-6
    pr(f'  {k}: pred(0,025) {Pa[k]:.10e} (silnik R4/pred - 1 {ra:+.2e}); pred(0,1) {Pb[k]:.10e} (silnik R2/pred - 1 {rb:+.2e}); '
       f'pred(0) {P0[k]:.10e}; {nap} apocentrow; start 0,025 vs 0,1 max |roznica| {dst:.1e}')
npa, npb, np0 = nu_of(Pa), nu_of(Pb), nu_of(P0)
okP &= abs(n4a - npa) <= 0.010e6 and abs(n2b - npb) <= 0.010e6 and conv <= 1e-6
pr(f'  nu_pred(0,025) = {npa/1e6:.4f} MHz (silnik {(n4a - npa)/1e6:+.4f}); nu_pred(0,1) = {npb/1e6:.4f} MHz (silnik {(n2b - npb)/1e6:+.4f}); '
   f'nu_pred(0) = {units(np0)}')
pr(f"H815p (GLOWNA 1: wezly <= 5e-6, nu <= 0,010 MHz, obie skale, zbieznosc <= 1e-6): {'TAK' if okP else 'NIE'} (zbieznosc max {conv:.1e})")
ratio = (nu_of(R2b) - n2a)/(npb - npa)
pr(f"H815e (skalowanie eps^2: [nu(0,1) - nu(0,025)] / [pred] w [0,99; 1,01]): {'TAK' if 0.99 <= ratio <= 1.01 else 'NIE'} ({ratio:.5f})")
dc = max(abs(PP[k]['lin']/A[(0.025, 2048, k)]['lin'] - 1) for k in NODES)
pr(f"H815c (p-Ps / o-Ps przy 0,025 i 2048: <= 1e-8): {'TAK' if dc <= 1e-8 else 'NIE'} (max {dc:.2e})")
nu0 = n4a - (n2b - n2a)/15
pr(f"H815z (granica eps -> 0: nu_0 = {units(nu0)}; pred(0) {np0/1e6:.4f}; roznica {(nu0 - np0)/1e6:+.4f} MHz; <= 0,010): "
   f"{'TAK' if abs(nu0 - np0) <= 0.010e6 else 'NIE'}")
# ---------------- H815o
c_step = n4a - nu_of({k: A[(0.025, 2048, k)]['lin'] for k in NODES}); c_eps = nu0 - n4a
nuP = nu_of({k: PP[k]['lin'] for k in NODES}) + c_step + c_eps; okO = True
for nm, nu in (('o-Ps (PHEN=2, eps -> 0)', nu0), ('p-Ps (PHEN=1, 2048 + poprawki kroku i eps z o-Ps)', nuP)):
    rho = nu/NU4; okO &= abs(rho - 1) <= 10*al**2
    pr(f'  {nm}: nu_orb = {units(nu)}; rho = {rho:.7f} = 1 {(rho - 1)/al**2:+.3f} alfa^2; delta_6 = {(nu - NU4)/1e6:+.4f} MHz '
       f'(pred(0) {(np0 - NU4)/1e6:+.4f}; jednocialowy +0,182); klasyczne 8/9 o {abs(rho - 8/9)/(10*al**2):.0f} progow')
pr(f"H815o (GLOWNA 2: |rho - 1| <= 10 alfa^2, o-Ps i p-Ps): {'TAK' if okO else 'NIE'}")
# ---------------- report: guard 0,05 r* on the 814 starts
G5, G0 = {}, {}
for k in NODES:
    tmp = os.path.join(D, 'out', f'_s814_{k}.err'); open(tmp, 'wb').write(tarbytes(T814, f'P5/F1_ops_{k}_s4096_5.err'))
    st = parse_start(tmp); os.remove(tmp); G5[k] = P15.omega(st, 0.05)[0]; G0[k] = P15.omega(st, 0.0)[0]
pr(f'  RAPORT oslona 0,05 r* (starty 814): nu_pred(0,05) - nu_pred(0) = {(nu_of(G5) - nu_of(G0))/1e6:+.4f} MHz '
   f'(814: silnik R4 - pred814 = -1,5208 MHz)')
allok = k1 and okV and okS and okP and 0.99 <= ratio <= 1.01 and dc <= 1e-8 and abs(nu0 - np0) <= 0.010e6 and okO
pr(f"DECYZJA: {'stala roznica 814 = zmiekczenie Plummera; <r^-3> trypletu Ps (n = 2, l = 1) w rzedzie alfa^4 = zmierzona roznica torusow orbity produkcyjnej (eps -> 0) -- regula Borna dla o-Ps przeklasyfikowana' if allok else 'regula Borna zostaje importem'}")

# Audit 814: Ps orbit under the PRODUCTION (conservative Coulomb + Darwin) force in the 699 measurement (CREM_PAIR_TORQUE_CONSERVATIVE=1):
# torus difference E(2,2) - E(2,1) from the apsidal precession against Born/811 (alpha^4) and against an independent integration of the
# same equations (pred814.py).  Usage: python3 ana814.py [runs dir]  (default out/runs) -> out/ana814.txt
import os, sys, math, tarfile
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
import pred814 as PRED
O = sys.argv[1] if len(sys.argv) > 1 else os.path.join(D, 'out', 'runs')
HEV, C = 4.135667696e-15, 299792458.0; al = 1/137.035999177; NU4 = 5474.642e6
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.4f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
def gauss(N):
    x, w = np.polynomial.legendre.leggauss(N); u = 0.75 + 0.25*x[::-1]; return 1/u, 0.25*w[::-1]
NODES = {N: dict(zip([f'u{N}{c}' for c in 'abcde'[:N]], zip(*gauss(N)))) for N in (3, 4, 5)}
def read(path, tag='mechB0'):   # = ana813.read (apse vertices as ana808)
    rows = []; has_sum = False; has_cons = False
    for ln in open(path):
        if ln.startswith(f'PAIR TORQUE799 {tag} sum'): has_sum = True
        if ln.startswith(f'PAIR CONS814 {tag} start'): has_cons = True
        if not ln.startswith(f'PAIR APSE808 {tag} '): continue
        w = ln.split(); tc = float(w[w.index('tc') + 1])
        r = [np.array([float(v) for v in w[w.index(f'r{q}') + 1:w.index(f'r{q}') + 4]]) for q in range(3)]
        rows.append((tc, r))
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
    lin = np.polyfit(T, P, 1)[0]; om0 = np.polyfit(T, P, 2)[1]
    Tk = np.diff(T); tm = 0.5*(T[1:] + T[:-1]); slope = np.polyfit(tm, np.log(Tk), 1)[0]
    return dict(om0=om0, lin=lin, slope=slope, n=len(T), nbad=nbad, has_sum=has_sum, has_cons=has_cons, path=path)
# ---------------- H814k
def tarbytes(tarpath, suffix):
    with tarfile.open(tarpath) as t:
        m = [x for x in t.getnames() if x.endswith(suffix)]; return t.extractfile(m[0]).read() if m else None
k1 = tarbytes(os.path.join(D, '..', 'section-813', 'runs', 'runs813.tar.xz'), 'P5/F1_ops_u5a_s1024_5.err') == open(f'{O}/K/F1_ops_u5a_s1024_5.err', 'rb').read()
k2 = tarbytes(os.path.join(D, '..', 'section-812', 'runs', 'runs812.tar.xz'), 'K/F1_proton_u3b_s2048_5.err') == open(f'{O}/K/F1_proton_u3b_s2048_5.err', 'rb').read()
pr(f"H814k (bez przelacznika: o-Ps u5a 1024 = 813, wodor u3b 2048 = 812, bajt w bajt): {'TAK' if k1 and k2 else 'NIE'} ({k1}, {k2})")
# ---------------- load
H = (1024, 2048, 4096)
OP = {(N, s, k): read(f'{O}/P{N}/F1_ops_{k}_s{s}_5.err') for N in (4, 5) for s in H for k in NODES[N]}
PP = {k: read(f'{O}/PP5/F0_pps_{k}_s2048_5.err') for k in NODES[5]}
HY = {F: read(f'{O}/H/F{F}_proton_u3b_s2048_5.err') for F in (1, 0)}
allf = list(OP.values()) + list(PP.values()) + list(HY.values())
okV = len(allf) == 34 and all(v['has_sum'] and v['has_cons'] and v['n'] >= 195 and v['nbad'] == 0 for v in allf)
pr(f"H814v (34 pliki: sum, CONS814, >= 195 apocentrow, 0 zlych wierzcholkow): {'TAK' if okV else 'NIE'} (plikow {len(allf)}, "
   f"min apocentrow {min(v['n'] for v in allf)}, zlych {sum(v['nbad'] for v in allf)})")
# ---------------- H814e
okE = True
for N in (4, 5):
    for k in NODES[N]:
        v = OP[(N, 4096, k)]; okE &= abs(v['slope']) <= 3e7
        pr(f'  o-Ps {k} 4096: nachylenie ln T {v["slope"]:+.3e} 1/s; Omega_0/Omega_lin - 1 {v["om0"]/v["lin"] - 1:+.2e}')
pr(f"H814e (|nachylenie ln T| <= 3e7 1/s, o-Ps 4096, 9 wezlow): {'TAK' if okE else 'NIE'}")
# ---------------- H814c
dc = max(abs(PP[k]['lin']/OP[(5, 2048, k)]['lin'] - 1) for k in NODES[5])
pr(f"H814c (p-Ps / o-Ps przy 2048, 5 wezlow, Omega_lin: <= 1e-8): {'TAK' if dc <= 1e-8 else 'NIE'} (max {dc:.2e})")
# ---------------- H814s, H814q
nu_of = lambda N, Om: sum(w*L**2*Om[k] for k, (L, w) in NODES[N].items())/(2*math.pi)
pobs = lambda y: math.log2((y[0] - y[1])/(y[1] - y[2])) if (y[0] - y[1])*(y[1] - y[2]) > 0 else float('nan')
okS = True; R = {}; R4n = {}
for N in (4, 5):
    R2 = {}; R4 = {}
    for k in NODES[N]:
        y = [OP[(N, s, k)]['lin'] for s in H]; p = pobs(y); okS &= 1.8 <= p <= 2.2
        R2[k] = y[2] + (y[2] - y[1])/3; r2p = y[1] + (y[1] - y[0])/3; R4[k] = R2[k] + (R2[k] - r2p)/15
        pr(f'  {N} wezly {k}: Omega_lin 1024 / 2048 / 4096 = ' + ' / '.join(f'{v:.10e}' for v in y) + f'; p = {p:.4f}; R4 {R4[k]:.10e}')
    ys = [nu_of(N, {k: OP[(N, s, k)]['lin'] for k in NODES[N]}) for s in H]; n2, n4 = nu_of(N, R2), nu_of(N, R4)
    okS &= abs(n4 - n2) <= 0.010e6; R[N] = (ys, n2, n4); R4n[N] = R4
    pr(f'  nu_orb {N} wezly 1024 / 2048 / 4096: ' + ' / '.join(f'{v/1e6:.4f}' for v in ys) + f' MHz; R2 {n2/1e6:.4f}, R4 {n4/1e6:.4f} MHz')
pr(f"H814s (p w [1,8; 2,2], |R4 - R2| <= 0,010 MHz, 4 i 5 wezlow): {'TAK' if okS else 'NIE'}")
dq = R[4][2] - R[5][2]
pr(f"H814q (|nu(4) - nu(5)| <= 0,010 MHz): {'TAK' if abs(dq) <= 0.010e6 else 'NIE'} ({dq/1e6:+.4f} MHz)")
# ---------------- H814n: independent integration
okN = True; OPRED = {}; conv = 0.0
for k in NODES[5]:
    sts = [PRED.parse_start(OP[(5, s, k)]['path']) for s in H]
    dst = max(max(np.max(np.abs(st[q] - sts[2][q])) for q in ('x1', 'x2', 'v1', 'v2')) for st in sts)
    a = PRED.omega(sts[2], rtol=1e-12)[0]; b, nap, Tr = PRED.omega(sts[2], rtol=1e-13); conv = max(conv, abs(b/a - 1)); OPRED[k] = b
    rr = R4n[5][k]/b - 1; okN &= abs(rr) <= 5e-6
    pr(f'  {k}: Omega_pred {b:.10e} rad/s ({nap} apocentrow, DOP853 1e-12/1e-13 {b/a - 1:+.1e}; start rozny miedzy krokami max {dst:.1e}); '
       f'silnik R4 / pred - 1 = {rr:+.2e}')
nup = nu_of(5, OPRED); okN &= abs(R[5][2] - nup) <= 0.010e6 and conv <= 1e-7
pr(f'  nu_pred (5 wezlow) = {units(nup)}; silnik R4 - pred = {(R[5][2] - nup)/1e6:+.4f} MHz')
pr(f"H814n (wezly <= 5e-6, nu <= 0,010 MHz, zbieznosc calki <= 1e-7): {'TAK' if okN else 'NIE'} (zbieznosc max {conv:.1e})")
# ---------------- H814h
okH = True
for F in (1, 0):
    with tarfile.open(os.path.join(D, '..', 'section-808', 'runs', 'runs808.tar.xz')) as t:
        m = [x for x in t.getnames() if x.endswith(f'O1/F{F}_proton_u3b_s2048_5.err')][0]
        tmp = os.path.join(D, 'out', f'_h808_F{F}.err'); open(tmp, 'wb').write(t.extractfile(m).read())
    ref = read(tmp); os.remove(tmp); rr = HY[F]['lin']/ref['lin'] - 1; okH &= abs(rr) <= 1e-6
    pr(f'  wodor F = {F} (PHEN={2 if F else 1}) u3b 2048: zachowawcza {HY[F]["lin"]:.10e}, opozniona 808 {ref["lin"]:.10e}, stosunek - 1 {rr:+.2e}')
pr(f"H814h (wodor zachowawcza vs opozniona, F = 1 i F = 0: <= 1e-6): {'TAK' if okH else 'NIE'}")
# ---------------- H814o
nuO = R[5][2]; c_o = R[5][2] - R[5][0][1]; nuP = nu_of(5, {k: PP[k]['lin'] for k in NODES[5]}) + c_o; okO = True
for nm, nu in (('o-Ps (PHEN=2, 5 wezlow, R4)', nuO), ('p-Ps (PHEN=1, 2048 + poprawka kroku o-Ps)', nuP)):
    rho = nu/NU4; okO &= abs(rho - 1) <= 10*al**2
    pr(f'  {nm}: nu_orb = {units(nu)}; rho = {rho:.7f} = 1 {(rho - 1)/al**2:+.3f} alfa^2; delta_6 = {(nu - NU4)/1e6:+.4f} MHz '
       f'(pred {(nup - NU4)/1e6:+.4f}; jednocialowy +0,182; 813 opozniona -1,506); klasyczne 8/9 o {abs(rho - 8/9)/(10*al**2):.0f} progow')
pr(f"H814o (GLOWNA: |rho - 1| <= 10 alfa^2, o-Ps i p-Ps): {'TAK' if okO else 'NIE'}")
allok = k1 and k2 and okV and okE and dc <= 1e-8 and okS and abs(dq) <= 0.010e6 and okN and okH and okO
pr(f"DECYZJA: {'<r^-3> trypletu Ps (n = 2, l = 1) w rzedzie alfa^4 = zmierzona roznica torusow orbity produkcyjnej -- regula Borna dla o-Ps przeklasyfikowana' if allok else 'regula Borna zostaje importem'}")

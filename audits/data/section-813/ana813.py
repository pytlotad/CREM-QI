# Audit 813: apsidal rate of the Ps orbit at the start torus (Omega_0, quadratic fit of the apse angle in time) instead of the run-average
# linear slope of 812; dissipation of the retarded mutual field (period drift) against -(4/3) nu K^3 <r^-4>/(mu^2 c^3).
# Inputs: out/runs (new, 5 nodes), ../section-812/runs/runs812.tar.xz, ../section-808/runs/runs808.tar.xz.   Usage: python3 ana813.py
import os, sys, math, tarfile, tempfile
import numpy as np
D = os.path.dirname(os.path.abspath(__file__))
NEW = sys.argv[1] if len(sys.argv) > 1 else os.path.join(D, 'out', 'runs')
TMP = tempfile.mkdtemp()
for sec, tarn in (('812', 'runs812.tar.xz'), ('808', 'runs808.tar.xz')):
    with tarfile.open(os.path.join(D, '..', f'section-{sec}', 'runs', tarn)) as t:
        t.extractall(os.path.join(TMP, sec), filter='data')
A812 = os.path.join(TMP, '812', 'runs'); A808 = os.path.join(TMP, '808', 'runs', 'O1')
HEV, C, h = 4.135667696e-15, 299792458.0, 6.62607015e-34
hb = h/(2*math.pi); al = 1/137.035999177; K = al*hb*C
me, mp = 9.1093837139e-31, 1.67262192595e-27
NU4 = 5474.642e6
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.4f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
def gauss(N):
    x, w = np.polynomial.legendre.leggauss(N); u = 0.75 + 0.25*x[::-1]; return 1/u, 0.25*w[::-1]
NODES = {N: dict(zip([f'u{N}{c}' for c in 'abcde'[:N]], zip(*gauss(N)))) for N in (3, 4, 5)}
def read(path, tag='mechB0'):
    rows = []; has_sum = False
    for ln in open(path):
        if ln.startswith(f'PAIR TORQUE799 {tag} sum'): has_sum = True
        if not ln.startswith(f'PAIR APSE808 {tag} '): continue
        w = ln.split(); tc = float(w[w.index('tc') + 1])
        r = [np.array([float(v) for v in w[w.index(f'r{q}') + 1:w.index(f'r{q}') + 4]]) for q in range(3)]
        rows.append((tc, r))
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
def pred_slope(m1, m2, lfrac):
    mu = m1*m2/(m1 + m2); nu = m1*m2/(m1 + m2)**2; a = 4*hb**2/(mu*K); e2 = 1 - lfrac**2
    r4 = (1 + e2/2)/(a**4*(1 - e2)**2.5); dEdt = -(4/3)*nu*K**3*r4/(mu**2*C**3)
    return 1.5*dEdt/(K/(2*a))
LF = lambda k, N: NODES[N][k][0]/2
# ---------------- load
O = {}   # (N, steps, node) -> dict, o-Ps
for N, sub in ((3, 'P3'), (4, 'P4')):
    for s in (1024, 2048, 4096):
        for k in NODES[N]: O[(N, s, k)] = read(f'{A812}/{sub}/F1_ops_{k}_s{s}_5.err')
for s in (1024, 2048, 4096):
    for k in NODES[5]: O[(5, s, k)] = read(f'{NEW}/P5/F1_ops_{k}_s{s}_5.err')
PP = {k: read(f'{NEW}/PP5/F0_pps_{k}_s2048_5.err') for k in NODES[5]}
nbad = sum(v['nbad'] for v in list(O.values()) + list(PP.values())); okfiles = all(v['has_sum'] and v['n'] >= 195 for v in list(O.values()) + list(PP.values()))
pr(f'Pliki Ps: {len(O) + len(PP)} (o-Ps {len(O)}, p-Ps {len(PP)}), wszystkie z linia sum i >= 195 apocentrami: {okfiles}; zlych wierzcholkow {nbad}')
# ---------------- H813d
okD = True
pr('H813d: nachylenie ln T [1/s], zmierzone / przewidziane -(4/3) nu K^3 <r^-4>/(mu^2 c^3) -> (3/2) dE/dt/|E|')
for N in (3, 4, 5):
    for k in NODES[N]:
        v = O[(N, 4096, k)]; p = pred_slope(me, me, LF(k, N)); r = v['slope']/p; okD &= 0.90 <= r <= 1.10
        pr(f'  o-Ps {k} 4096 (LFRAC {LF(k, N):.5f}): {v["slope"]:+.4e} / {p:+.4e} = {r:.4f}')
for k in NODES[5]:
    v = PP[k]; p = pred_slope(me, me, LF(k, 5)); pr(f'  p-Ps {k} 2048 (raport): {v["slope"]:+.4e} / {p:+.4e} = {v["slope"]/p:.4f}')
HYD = []
for F in (1, 0):
    for k in NODES[3]: HYD.append((f'808 O1 F = {F} {k}', read(f'{A808}/F{F}_proton_{k}_s2048_5.err'), LF(k, 3)))
HYD.append(('812 K F = 1 u3b', read(f'{A812}/K/F1_proton_u3b_s2048_5.err'), LF('u3b', 3)))
for nm, v, lf in HYD:
    pPs = pred_slope(me, me, lf); pH = pred_slope(mp, me, lf); okD &= abs(v['slope']) <= 0.01*abs(pPs)
    pr(f'  wodor {nm}: {v["slope"]:+.4e} (przewidywanie wodoru {pH:+.4e}, Ps przy tym L {pPs:+.4e}; |zm|/|Ps| {abs(v["slope"]/pPs):.2e})')
pr(f"H813d (Ps 4096: zmierzone/przewidziane w [0,90; 1,10]; wodor |zm| <= 0,01 Ps): {'TAK' if okD else 'NIE'}")
# ---------------- H813h
dh = max(abs(v['om0']/v['lin'] - 1) for nm, v, lf in HYD[:6])
for nm, v, lf in HYD: pr(f'  wodor {nm}: Omega_0 {v["om0"]:.10e}, Omega_lin {v["lin"]:.10e}, stosunek - 1 {v["om0"]/v["lin"] - 1:+.2e}')
pr(f"H813h (wodor 808 O1, F = 1 i F = 0: |Omega_0/Omega_lin - 1| <= 1e-6): {'TAK' if dh <= 1e-6 else 'NIE'} (max {dh:.2e})")
# ---------------- H813s, H813q
nu_of = lambda N, Om: sum(w*L**2*Om[k] for k, (L, w) in NODES[N].items())/(2*math.pi)
pobs = lambda y: math.log2((y[0] - y[1])/(y[1] - y[2])) if (y[0] - y[1])*(y[1] - y[2]) > 0 else float('nan')
H = (1024, 2048, 4096); okS = True; R = {}
for N in (3, 4, 5):
    R2 = {}; R4 = {}
    for k in NODES[N]:
        y = [O[(N, s, k)]['om0'] for s in H]; p = pobs(y); okS &= 1.8 <= p <= 2.2
        R2[k] = y[2] + (y[2] - y[1])/3; r2p = y[1] + (y[1] - y[0])/3; R4[k] = R2[k] + (R2[k] - r2p)/15
        pr(f'  {N} wezly {k}: Omega_0 1024 / 2048 / 4096 = ' + ' / '.join(f'{v:.10e}' for v in y) + f'; p = {p:.4f}; (Omega_lin 4096 {O[(N, 4096, k)]["lin"]:.10e})')
    ys = [nu_of(N, {k: O[(N, s, k)]['om0'] for k in NODES[N]}) for s in H]; n2, n4 = nu_of(N, R2), nu_of(N, R4)
    okS &= abs(n4 - n2) <= 0.01e6; R[N] = (ys, n2, n4)
    pr(f'  nu_orb(Omega_0) {N} wezly 1024 / 2048 / 4096: ' + ' / '.join(f'{v/1e6:.4f}' for v in ys) + f' MHz; R2 {n2/1e6:.4f}, R4 {n4/1e6:.4f} MHz')
pr(f"H813s (p w [1,8; 2,2] w kazdym wezle, |R4 - R2| <= 0,01 MHz dla 3, 4, 5 wezlow): {'TAK' if okS else 'NIE'}")
d35 = R[3][2] - R[5][2]; d45 = R[4][2] - R[5][2]; okQ = abs(d35) <= 0.010e6 and abs(d45) <= 0.010e6
pr(f"H813q (|nu(3) - nu(5)|, |nu(4) - nu(5)| <= 0,010 MHz, Omega_0, R4): {'TAK' if okQ else 'NIE'} ({d35/1e6:+.4f}, {d45/1e6:+.4f} MHz)")
# ---------------- H813c
dc = max(abs(PP[k]['om0']/O[(5, 2048, k)]['om0'] - 1) for k in NODES[5])
pr(f"H813c (p-Ps / o-Ps przy 2048, 5 wezlow, Omega_0: <= 1e-8): {'TAK' if dc <= 1e-8 else 'NIE'} (max {dc:.2e})")
# ---------------- H813o
nuO = R[5][2]; c_o = R[5][2] - R[5][0][1]; nuP = nu_of(5, {k: PP[k]['om0'] for k in NODES[5]}) + c_o
a2 = al**2; okO = True
for nm, nu in (('o-Ps (PHEN=2, 5 wezlow, R4, Omega_0)', nuO), ('p-Ps (PHEN=1, 2048 + poprawka kroku o-Ps)', nuP)):
    rho = nu/NU4; okO &= abs(rho - 1) <= 10*a2
    pr(f'  {nm}: nu_orb = {units(nu)}; rho_0 = {rho:.7f} = 1 + {(rho - 1)/a2:.3f} alfa^2; delta_6 = {(nu - NU4)/1e6:+.4f} MHz '
       f'(jednocialowy +0,182; 812 liniowo +1,278); klasyczne 8/9 odlegle o {abs(rho - 8/9)/(10*a2):.0f} progow')
pr(f"H813o (GLOWNA: |rho_0 - 1| <= 10 alfa^2 = {10*a2:.2e}, o-Ps i p-Ps): {'TAK' if okO else 'NIE'}")
nl = {N: nu_of(N, {k: O[(N, 4096, k)]['lin'] for k in NODES[N]}) for N in (3, 4, 5)}
pr('  raport: nu_orb z Omega_lin przy 4096 (jak 812): ' + ', '.join(f'{N} wezly {nl[N]/1e6:.4f}' for N in (3, 4, 5)) + ' MHz')
allok = okD and dh <= 1e-6 and okS and okQ and dc <= 1e-8 and okO
pr(f"DECYZJA: {'<r^-3> trypletu Ps (n = 2, l = 1) w rzedzie alfa^4 = zmierzona roznica torusow orbity zachowawczej -- przeklasyfikowanie w Model.md' if allok else 'regula Borna zostaje importem'}")

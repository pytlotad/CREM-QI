# Audit 809: (1) quadrature control of 808 at the step limit -- 4 Gauss nodes in u = 1/L at 1024 / 2048 (808 O3) / 4096 against the
# 3-node limit of 808; (2) C image of the antihydrogen anomaly runs (hydrogen with the heavy spin reversed, CREM_PAIR_TORQUE_HEAVY_FLIP=1).
# usage: python3 ana809.py [809 runs dir] [808 runs dir]   (808: runs808.tar.xz of section-808 unpacked; default ../section-808/out/runs)
import os, sys, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
from per_orbit import om
O = sys.argv[1] if len(sys.argv) > 1 else os.path.join(D, 'out', 'runs')
P = sys.argv[2] if len(sys.argv) > 2 else os.path.join(D, '..', 'section-808', 'out', 'runs')
HEV, C = 4.135667696e-15, 299792458.0
NU_NAT = 10969.045e6; S_REF = math.hypot(0.015e6, 0.020e6)
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.4f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm = {f/NU_NAT:.6f} pomiaru'
def gauss(N):
    x, w = np.polynomial.legendre.leggauss(N); u = 0.75 + 0.25*x[::-1]; return 1/u, 0.25*w[::-1]   # L ascending
N3 = dict(zip(('u3a', 'u3b', 'u3c'), zip(*gauss(3)))); N4 = dict(zip(('u4a', 'u4b', 'u4c', 'u4d'), zip(*gauss(4))))
def apse(path, tag='mechB0'):
    """Omega_aps [rad/s] from the parabola vertex r(theta) at every apoapsis (as ana808.py); (Omega, n_apo, n_bad_vertex, has_sum)."""
    rows = []; has_sum = False
    for ln in open(path):
        if ln.startswith(f'PAIR TORQUE799 {tag} sum'): has_sum = True
        if not ln.startswith(f'PAIR APSE808 {tag} '): continue
        w = ln.split(); tc = float(w[w.index('tc') + 1])
        r = [np.array([float(v) for v in w[w.index(f'r{q}') + 1:w.index(f'r{q}') + 4]]) for q in range(3)]
        rows.append((tc, r))
    if not rows: return float('nan'), 0, 0, has_sum
    r0 = rows[0][1]; Lh = np.cross(r0[0], r0[2]); Lh /= np.linalg.norm(Lh)
    e1 = r0[1]/np.linalg.norm(r0[1]); e1 -= Lh*np.dot(Lh, e1); e1 /= np.linalg.norm(e1); e2 = np.cross(Lh, e1)
    T = []; Ph = []; nbad = 0; prev = None
    for tc, r in rows:
        th = np.unwrap([math.atan2(np.dot(x, e2), np.dot(x, e1)) for x in r]); rho = [np.linalg.norm(x) for x in r]
        x = th - th[1]; A = np.vstack([np.ones(3), x, x**2]).T; a0, b, cq = np.linalg.solve(A, rho)
        xs = -b/(2*cq); h = max(abs(x[0]), abs(x[2]))
        if not (cq < 0 and x[0] - 0.5*h <= xs <= x[2] + 0.5*h): nbad += 1
        phi = th[1] + xs
        if prev is not None: phi += 2*math.pi*round((prev - phi)/(2*math.pi))
        prev = phi; T.append(tc); Ph.append(phi)
    k, b0 = np.polyfit(np.array(T), np.array(Ph), 1)
    return k, len(T), nbad, has_sum
f = lambda base, sub, F, t, node, h, s=5: f'{base}/{sub}/F{F}_{t}_{node}_s{h}_{s}.err'
okV = True; nnew = 0; minapo = 10**9; nbadv = 0; rejC = 0
def AP(path, new=True):
    global okV, nnew, minapo, nbadv
    k, n, nb, hs = apse(path)
    if new: nnew += 1; minapo = min(minapo, n); nbadv += nb; okV &= hs and n >= 195 and nb == 0
    return k
def OM(path, new=True):
    global okV, nnew, minapo, rejC
    hs = any(ln.startswith('PAIR TORQUE799 mechB0 sum') for ln in open(path)); o, r, _, nT = om(path)
    if new: nnew += 1; rejC += r; minapo = min(minapo, nT); okV &= hs and r == 0 and nT >= 195
    return o
nu_of = lambda nodes, Om: sum(w*L**2*Om[k] for k, (L, w) in nodes.items())/(2*math.pi)
# --- 4 nodes: 1024, 4096 (Q1, F = 1 hydrogen), 2048 (808 O3, F = 1 hydrogen), 2048 F = 0 (Q2)
y4 = {k: [AP(f(O, 'Q1', 1, 'proton', k, 1024)), AP(f(P, 'O3', 1, 'proton', k, 2048), False), AP(f(O, 'Q1', 1, 'proton', k, 4096))] for k in N4}
q2 = {k: AP(f(O, 'Q2', 0, 'proton', k, 2048)) for k in N4}
# --- 3 nodes from 808 (O2 1024, O1 2048, O5 4096; F = 1 hydrogen)
y3 = {k: [AP(f(P, 'O2', 1, 'proton', k, 1024), False), AP(f(P, 'O1', 1, 'proton', k, 2048), False), AP(f(P, 'O5', 1, 'proton', k, 4096), False)] for k in N3}
# --- C image: hydrogen flipped (C1) vs antihydrogen 808 A1 (1 - F)
Cf = {(F, g): OM(f(O, 'C1', F, 'proton', 'flip' + g, 2048)) for F in (1, 0) for g in ('gE', 'g2')}
Ah = {(F, g): OM(f(P, 'A1', F, 'antiproton', g, 2048), False) for F in (1, 0) for g in ('gE', 'g2')}
Hh = {(F, g): OM(f(P, 'A1', F, 'proton', g, 2048), False) for F in (1, 0) for g in ('gE', 'g2')}
okV &= nnew == 16
pr(f"H809v (16 nowych plikow: linia mechB0 sum; Q >= 195 apocentrow APSE808, 0 zlych wierzcholkow; C 0 odrzuconych, >= 195 apocentrow): "
   f"{'TAK' if okV else 'NIE'} (nowych plikow {nnew}, apocentrow min {minapo}, zlych wierzcholkow {nbadv}, odrzuconych C {rejC})")
def lim(y):
    p = math.log2((y[0] - y[1])/(y[1] - y[2])) if (y[0] - y[1])*(y[1] - y[2]) > 0 else float('nan')
    R2 = y[2] + (y[2] - y[1])/3; R2p = y[1] + (y[1] - y[0])/3; return p, R2, R2 + (R2 - R2p)/15
okS = True; R2_4 = {}; R4_4 = {}; R4_3 = {}; R2_3 = {}
for k, (L, w) in N4.items():
    p, R2_4[k], R4_4[k] = lim(y4[k]); okS &= 1.9 <= p <= 2.1
    pr(f'  {k} L = {L:.6f}: Omega_aps 1024 / 2048 / 4096 = ' + ' / '.join(f'{v:.10e}' for v in y4[k]) + f'; p = {p:.4f}; R2 {R2_4[k]:.10e}, R4 {R4_4[k]:.10e} rad/s')
for k in N3: _, R2_3[k], R4_3[k] = lim(y3[k])
n4R4, n4R2 = nu_of(N4, R4_4), nu_of(N4, R2_4); n3R4, n3R2 = nu_of(N3, R4_3), nu_of(N3, R2_3)
n4h = {h: nu_of(N4, {k: y4[k][i] for k in N4}) for i, h in enumerate((1024, 2048, 4096))}
n3h = {h: nu_of(N3, {k: y3[k][i] for k in N3}) for i, h in enumerate((1024, 2048, 4096))}
okS &= abs(n4R4 - n4R2) <= 0.02e6
pr('  nu_orb 4 wezly 1024 / 2048 / 4096: ' + ' / '.join(f'{n4h[h]/1e6:.4f}' for h in n4h) + f' MHz; R2 {n4R2/1e6:.4f}, R4 {n4R4/1e6:.4f}')
pr('  nu_orb 3 wezly 1024 / 2048 / 4096: ' + ' / '.join(f'{n3h[h]/1e6:.4f}' for h in n3h) + f' MHz; R2 {n3R2/1e6:.4f}, R4 {n3R4/1e6:.4f} (808)')
pr('  roznica 3 - 4 wezly: 1024 / 2048 / 4096 / R2 / R4 = ' + ' / '.join(f'{(n3h[h] - n4h[h])/1e6:+.4f}' for h in n3h)
   + f' / {(n3R2 - n4R2)/1e6:+.4f} / {(n3R4 - n4R4)/1e6:+.4f} MHz')
pr(f"H809s (p w [1,9; 2,1] w kazdym z 4 wezlow, |nu4(R4) - nu4(R2)| <= 0,02 MHz): {'TAK' if okS else 'NIE'} (|R4 - R2| = {abs(n4R4 - n4R2)/1e6:.4f} MHz)")
dc = max(abs(q2[k]/y4[k][1] - 1) for k in N4)
okC = dc <= 1e-9
pr(f"H809c (Q2 F = 0 wobec 808 O3 F = 1, 2048, <= 1e-9 wzgl.): {'TAK' if okC else 'NIE'} (max {dc:.2e})")
dq = n3R4 - n4R4; okQ = abs(dq) <= 0.010e6
pr(f"H809q (GLOWNA dla 808; |nu3(R4) - nu4(R4)| <= 0,010 MHz): {'TAK' if okQ else 'NIE'} (nu3 {n3R4/1e6:.4f}, nu4 {n4R4/1e6:.4f}, roznica {dq/1e6:+.4f} MHz)")
pr('Symetria C: wodor z odwroconym spinem ciezkiej (C1, F) wobec antywodoru 808 A1 (1 - F):')
dC = 0.0
for F in (1, 0):
    for g in ('gE', 'g2'):
        d = Cf[(F, g)]/Ah[(1 - F, g)] - 1; dC = max(dC, abs(d))
        pr(f'  F = {F} {g}: wodor odwr. {Cf[(F, g)]:.10e} Hz, antywodor F = {1 - F} {Ah[(1 - F, g)]:.10e} Hz, wzgl. {d:+.2e}; wodor bez odwr. {Hh[(F, g)]:.10e}')
nAf = {F: 1.125*(Cf[(F, 'gE')] - Cf[(F, 'g2')]) for F in (1, 0)}
nAa = {F: 1.125*(Ah[(F, 'gE')] - Ah[(F, 'g2')]) for F in (1, 0)}
nAh = {F: 1.125*(Hh[(F, 'gE')] - Hh[(F, 'g2')]) for F in (1, 0)}
dA = max(abs(nAf[F] - nAa[1 - F]) for F in (1, 0))
for F in (1, 0):
    pr(f'  nu_A F = {F}: wodor odwr. {nAf[F]/1e6:.5f}, antywodor (F = {1 - F}) {nAa[1 - F]/1e6:.5f}, wodor {nAh[F]/1e6:.5f} MHz; '
       f'srednia po s1 = +-u (wodor, wodor odwr.) {0.5*(nAh[F] + nAf[F])/1e6:.5f} MHz')
okCC = dC <= 1e-9 and dA <= 1e-4*1e6
pr(f"H809C (<= 1e-9 wzgl. i nu_A <= 1e-4 MHz): {'TAK' if okCC else 'NIE'} (max wzgl. {dC:.2e}, max |d nu_A| {dA/1e6:.2e} MHz)")
# diagnosis: X = Omega(flipped) - Omega(hydrogen) vs 2 Omega_L of the heavy dipole field (Kepler <r^-3>, a_0 with m_r)
mu0 = 4e-7*math.pi*1.00000000055; e = 1.602176634e-19; me = 9.1093837139e-31; mp = 1.67262192595e-27; mup = 1.41060679545e-26
a0 = 5.29177210544e-11*(1 + me/mp); cth = 0.7696371619001129
OmL = (e/me)*mu0/(4*math.pi)/(8*1.5**3*a0**3)*mup*cth
for F in (1, 0):
    X = Cf[(F, 'g2')] - Hh[(F, 'g2')]
    pr(f'  DIAGNOZA F = {F}: X = {X:.6e} Hz = {2*math.pi*X:.6e} rad/s; 2 Omega_L = {2*OmL:.6e} rad/s; iloraz {2*math.pi*X/(2*OmL):.5f}; '
       f'X(g_e)/X(2) = {(Cf[(F, "gE")] - Hh[(F, "gE")])/X:.7f}')
s_o = max(abs(n4R4 - n4R2), abs(dq), 0.005e6); s_A = 0.010e6
fin = {}
for F in (1, 0):
    na = 0.5*(nAh[F] + nAa[F])
    nu = n4R4 + na; sig = math.sqrt(s_o**2 + s_A**2 + S_REF**2); z = (nu - NU_NAT)/sig
    pr(f'=== F = {F} (PHEN={F + 1}) ===')
    pr(f'  nu_orb = nu4(R4) = {units(n4R4)} (sigma_o {s_o/1e6:.4f} MHz); orbita kanalowo niezalezna (H808c, H809c)')
    pr(f'  nu_A = {na/1e6:.4f} MHz (808, para wodor + antywodor)')
    pr(f'  nu = {units(nu)}; nu - nu_nat = {(nu - NU_NAT)/1e6:+.4f} MHz = {z:+.2f} sigma (sigma {sig/1e6:.4f} MHz)')
    pr(f"  H809n F = {F} (|nu - nu_nat| <= 3 sigma): {'TAK' if abs(z) <= 3 else 'NIE'}")
    fin[F] = abs(z) <= 3; r = nu/NU_NAT
    pr(f'  o-Ps przez czynnik z wodoru: {r*0.9998:.5f}-{r*0.9999:.5f} multipletu (Hagena i in., PRL 71, 2887 (1993)); p-Ps bez czesci spinowej')
ok808 = okQ and okV and okS and okC
pr(f'WYNIK: 808 "orbita niesie 5/8 alfa^2" {"POTWIERDZONY" if ok808 else "NADAL NIEPOTWIERDZONY"}; H809C {"TAK" if okCC else "NIE"}; '
   f'H809n F = 1 {"TAK" if fin[1] else "NIE"}, F = 0 {"TAK" if fin[0] else "NIE"}')
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)).')
pr('Import QED 0,000 +- 0,020 MHz: Blumer i in., arXiv:2412.19580, Tab. 1; Horbatsch i Hessels, PRA 93, 022513 (2016).')

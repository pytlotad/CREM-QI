# Audit 808: hydrogen 2P3/2 - 2P1/2 from the apsidal precession of the mechanical orbit integrated over L along J_r + L = 2 hbar
# (Gauss-Legendre in u = 1/L, Sommerfeld k = j + 1/2) + the anomalous part of the BMT spin precession (g_e vs g = 2 at L = 1.5).
# usage: python3 ana808.py [runs dir]   (default out/runs; O1, O2, O3, O4, A1, A2 as in runs/all.sh; O5 (4096) from runs/o5.sh, POPRAWKA 2)
import os, sys, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
from per_orbit import om
O = sys.argv[1] if len(sys.argv) > 1 else os.path.join(D, 'out', 'runs')
HEV, C = 4.135667696e-15, 299792458.0
NU_NAT = 10969.045e6; S_REF = math.hypot(0.015e6, 0.020e6)
LEAD, EM = 10943.3239e6, 10943.6881e6        # Pauli with m_r (own calculation), Dirac E_M with m_r (Blumer et al., Tab. 1)
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.4f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm = {f/NU_NAT:.6f} pomiaru'
def gauss(N):
    x, w = np.polynomial.legendre.leggauss(N); u = 0.75 + 0.25*x[::-1]; return 1/u, 0.25*w[::-1]   # L ascending
N3 = dict(zip(('u3a', 'u3b', 'u3c'), zip(*gauss(3)))); N4 = dict(zip(('u4a', 'u4b', 'u4c', 'u4d'), zip(*gauss(4))))
TYPES = ('proton', 'antiproton')
bad = []
def apse(path, tag='mechB0'):
    """Omega_aps [rad/s] from the parabola vertex r(theta) at every apoapsis; returns (Omega, n_apo, n_bad_vertex, rms residual rad, has_sum)."""
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
        xs = -b/(2*cq); h = max(abs(x[0]), abs(x[2]))
        if not (cq < 0 and x[0] - 0.5*h <= xs <= x[2] + 0.5*h): nbad += 1
        phi = th[1] + xs
        if prev is not None: phi += 2*math.pi*round((prev - phi)/(2*math.pi))   # nearest branch to the previous apse (slow drift)
        prev = phi; T.append(tc); P.append(phi)
    T = np.array(T); P = np.array(P); k, b0 = np.polyfit(T, P, 1); res = P - (k*T + b0)
    return k, len(T), nbad, float(np.sqrt(np.mean(res**2))), has_sum
def f(sub, F, t, node, h, s=5): return f'{O}/{sub}/F{F}_{t}_{node}_s{h}_{s}.err'
okV = True; nfiles = 0; minapo = 10**9; rejA = 0; nbadv = 0
def AP(path):
    global okV, nfiles, minapo, nbadv
    k, n, nb, rms, hs = apse(path); nfiles += 1; minapo = min(minapo, n); nbadv += nb; okV &= hs and n >= 195 and nb == 0; return k, rms
def OM(path):
    global okV, nfiles, rejA, minapo
    hs = any(ln.startswith('PAIR TORQUE799 mechB0 sum') for ln in open(path)); o, r, _, nT = om(path); nfiles += 1
    rejA += r; minapo = min(minapo, nT); okV &= hs and r == 0 and nT >= 195; return o
nu_of = lambda nodes, Om: sum(w*L**2*Om[k] for k, (L, w) in nodes.items())/(2*math.pi)
# O1
Om1 = {(F, t, k): AP(f('O1', F, t, k, 2048)) for F in (1, 0) for t in TYPES for k in N3}
Om3 = {k: AP(f('O3', 1, 'proton', k, 2048))[0] for k in N4}
Om4 = AP(f('O4', 1, 'proton', 'u3b', 2048, 6))[0]
Om2 = {(h, k): AP(f('O2', 1, 'proton', k, h))[0] for h in (512, 1024) for k in N3}
for k in N3: Om2[(2048, k)] = Om1[(1, 'proton', k)][0]; Om2[(4096, k)] = AP(f('O5', 1, 'proton', k, 4096))[0]
nuA = {}; A2 = {}
for F in (1, 0):
    for t in TYPES: nuA[(F, t)] = 1.125*(OM(f('A1', F, t, 'gE', 2048)) - OM(f('A1', F, t, 'g2', 2048)))
A1024 = 1.125*(OM(f('A2', 1, 'proton', 'gE', 1024)) - OM(f('A2', 1, 'proton', 'g2', 1024)))
okV &= nfiles == 36
pr(f"H808v (36 plikow: linia mechB0 sum; >= 195 apocentrow; wierzcholki parabol w przedziale; 0 odrzuconych w A): {'TAK' if okV else 'NIE'} "
   f"(plikow {nfiles}, apocentrow min {minapo}, zlych wierzcholkow {nbadv}, odrzuconych A {rejA})")
pr('Omega_aps przy 2048 krokach [rad/s] i rms reszt dopasowania [rad]:')
for k, (L, w) in N3.items():
    pr(f'  {k} L = {L:.6f}: ' + '; '.join(f'F{F} {t} {Om1[(F, t, k)][0]:.10e} (rms {Om1[(F, t, k)][1]:.1e})' for F in (1, 0) for t in TYPES))
dev = max(abs(Om1[(F, t, k)][0]/Om1[(1, 'proton', k)][0] - 1) for F in (1, 0) for t in TYPES for k in N3)
ds = abs(Om4/Om1[(1, 'proton', 'u3b')][0] - 1)
pr(f"H808c (kanaly <= 1e-9, ziarno 6 <= 1e-7): {'TAK' if dev <= 1e-9 and ds <= 1e-7 else 'NIE'} (max kanaly {dev:.2e}, ziarno 6 {ds:.2e})")
H = (512, 1024, 2048, 4096); okS = True
pobs = lambda y: math.log2((y[0] - y[1])/(y[1] - y[2])) if (y[0] - y[1])*(y[1] - y[2]) > 0 else float('nan')
R2 = {}; R4 = {}
for k in N3:
    y = [Om2[(h, k)] for h in H]; pa, pb = pobs(y[0:3]), pobs(y[1:4]); okS &= 1.9 <= pa <= 2.1 and 1.9 <= pb <= 2.1
    R2[k] = y[3] + (y[3] - y[2])/3; r2p = y[2] + (y[2] - y[1])/3; R4[k] = R2[k] + (R2[k] - r2p)/15
    pr(f'  krok {k}: Omega 512 / 1024 / 2048 / 4096 = ' + ' / '.join(f'{v:.10e}' for v in y) + f'; p = {pa:.4f} / {pb:.4f}; R2 {R2[k]:.10e}, R4 {R4[k]:.10e}')
ys = {h: nu_of(N3, {k: Om2[(h, k)] for k in N3}) for h in H}
nR2, nR4 = nu_of(N3, R2), nu_of(N3, R4); c_o = nR4 - ys[2048]; s_o = max(abs(nR4 - nR2), 0.005e6); okS &= abs(nR4 - nR2) <= 0.02e6
pr('  nu_orb (F = 1 wodor) 512 / 1024 / 2048 / 4096: ' + ' / '.join(f'{ys[h]/1e6:.4f}' for h in H) + f' MHz; R2 {nR2/1e6:.4f}, R4 {nR4/1e6:.4f} MHz; c_o = {c_o/1e6:+.4f} MHz, sigma_o {s_o/1e6:.4f}')
pr(f"H808s (p w [1,9; 2,1] w kazdym wezle z obu trojek, |R4 - R2| <= 0,02 MHz): {'TAK' if okS else 'NIE'}")
q3 = ys[2048]; q4 = nu_of(N4, Om3)
pr(f"H808q (|nu_orb(3) - nu_orb(4)| <= 0,010 MHz, F = 1 wodor, 2048): {'TAK' if abs(q3 - q4) <= 0.010e6 else 'NIE'} (3 wezly {q3/1e6:.4f}, 4 wezly {q4/1e6:.4f}, roznica {(q3 - q4)/1e6:+.4f} MHz)")
c_A = nuA[(1, 'proton')] - A1024; s_A = max(abs(c_A), 0.010e6)
pr(f'  nu_A (F = 1 wodor) 1024 {A1024/1e6:.4f}, 2048 {nuA[(1, "proton")]/1e6:.4f} MHz; c_A = {c_A/1e6:+.4f}, sigma_A {s_A/1e6:.4f} MHz')
fin = {}
for F in (1, 0):
    pr(f'=== F = {F} (PHEN={F + 1}); pary (wodor + antywodor)/2 ===')
    no = nu_of(N3, {k: 0.5*(Om1[(F, 'proton', k)][0] + Om1[(F, 'antiproton', k)][0]) for k in N3}) + c_o
    na = 0.5*(nuA[(F, 'proton')] + nuA[(F, 'antiproton')])
    fr = (no - LEAD)/(EM - LEAD); nu = no + na; sig = math.sqrt(s_o**2 + s_A**2 + S_REF**2); z = (nu - NU_NAT)/sig
    pr(f'  nu_orb (granica kroku) = {units(no)}; wodor {nu_of(N3, {k: Om1[(F, "proton", k)][0] for k in N3})/1e6:.4f}, antywodor {nu_of(N3, {k: Om1[(F, "antiproton", k)][0] for k in N3})/1e6:.4f} MHz przy 2048')
    pr(f'  nu_A = {na/1e6:.4f} MHz (wodor {nuA[(F, "proton")]/1e6:.4f}, antywodor {nuA[(F, "antiproton")]/1e6:.4f})')
    pr(f"  H808o F = {F} (|f - 1| <= 0,1): {'TAK' if abs(fr - 1) <= 0.1 else 'NIE'} (f = {fr:.4f}; sigma_f = {s_o/(EM - LEAD):.3f})")
    pr(f'  nu = nu_orb + nu_A + delta_QED(0) = {units(nu)}; nu - nu_nat = {(nu - NU_NAT)/1e6:+.4f} MHz = {z:+.2f} sigma (sigma {sig/1e6:.4f} MHz)')
    pr(f"  H808n F = {F} (|nu - nu_nat| <= 3 sigma): {'TAK' if abs(z) <= 3 else 'NIE'}")
    fin[F] = abs(z) <= 3; r = nu/NU_NAT
    pr(f'  o-Ps przez czynnik z wodoru: {r*0.9998:.5f}-{r*0.9999:.5f} multipletu (Hagena i in., PRL 71, 2887 (1993)); p-Ps bez czesci spinowej')
allok = fin[1] and fin[0] and okV and okS and dev <= 1e-9 and ds <= 1e-7 and abs(q3 - q4) <= 0.010e6
pr(f'WYNIK: {"orbita niesie 5/8 alfa^2 (H808n TAK dla F = 0 i F = 1, kontrole TAK)" if allok else "NIE potwierdzone (patrz hipotezy)"}')
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)).')
pr('Import QED 0,000 +- 0,020 MHz: Blumer i in., arXiv:2412.19580, Tab. 1; Horbatsch i Hessels, PRA 93, 022513 (2016).')

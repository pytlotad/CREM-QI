# Audit 761 (offline): fine structure from spin precession (AC + Thomas) on the EBK torus; F basin under classical precession;
# spin phase per radial period vs pi.  -> out/ana761.txt   (see prereg761.txt)
import os, sys, math, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(D, '..', 'section-733'))
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()): import ana733 as w733          # EBK wave of principle I (u_ebk2), exact u_qm, grid
out = []; pr = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))
ALPHA = 7.2973525643e-3; AE = 1.15965218059e-3; MP_ME = 1836.15267343; MU = 1 / (1 + 1 / MP_ME)
HARTREE_HZ = 6.579683920502e15; HARTREE_EV = 27.211386245988; c_cm = 2.99792458e10
MEAS = 9911.200 + 1057.845                           # MHz (Hagley-Pipkin 1994 + Lundeen-Pipkin 1981)
def units(mhz): return f"{mhz:.3f} MHz = {mhz*1e6*HARTREE_EV/HARTREE_HZ*1e6:.4f} µeV = {c_cm/(mhz*1e6):.4f} cm"
n, l = 2, 1; L = l + 0.5
# --- <r^-3> variants (pair a.u., mu = 1) ---
a = n * n; e = math.sqrt(1 - L * L / (n * n)); u = np.linspace(0, 2 * np.pi, 2000001)
r = a * (1 - e * np.cos(u)); w = (1 - e * np.cos(u))                         # dt ~ (1 - e cos u) du
r3_orbit = np.trapezoid(w / r**3, u) / np.trapezoid(w, u)
rg = w733.grid(n); ue, ce = w733.u_ebk2(n, l, rg); uq = w733.u_qm(n, l, rg)
r3_ebk = np.trapezoid(ue * ue / rg**3, rg); r3_qm = np.trapezoid(uq * uq / rg**3, rg)
pr(f"<r^-3> (n 2, l 1): orbita EBK {r3_orbit:.8f} (1/(n^3 L^3) = {1/(n**3*L**3):.8f}); fala EBK 733 {r3_ebk:.8f} (warunek fazy {ce:.6f}); "
   f"QM {r3_qm:.8f} (1/24 = {1/24:.8f})")
pr(f"  stosunki do QM: orbita {r3_orbit/r3_qm:.6f}, fala EBK {r3_ebk/r3_qm:.6f}")
# --- part A ---
dLS_ebk = ((2.0)**2 - (1.0)**2) / 2; dLS_qm = (3.75 - 0.75) / 2
pr(f"roznica L.S j 3/2 - 1/2: EBK (|J| 2 i 1) {dLS_ebk}, QM {dLS_qm}")
def split(r3, g, recoil=True):
    m_M = 1 / MP_ME
    fac = (MU**3) * ((g - 1) + (m_M * g if recoil else 0)) if recoil else (g - 1)
    return 0.5 * ALPHA**2 * r3 * dLS_ebk * fac * HARTREE_HZ / 1e6      # MHz
g0, g1 = 2.0, 2 * (1 + AE)
V = [('V1 orbita EBK, g = 2', split(r3_orbit, g0)), ('V2 orbita EBK, g z a_e', split(r3_orbit, g1)),
     ('V3 fala EBK 733, g z a_e', split(r3_ebk, g1)), ('V4 QM <r^-3>, g z a_e', split(r3_qm, g1)), ('V4b QM <r^-3>, g = 2', split(r3_qm, g0))]
def somm(nr, k, mu):
    return mu / math.sqrt(1 + (ALPHA / (nr + math.sqrt(k * k - ALPHA * ALPHA)))**2) / ALPHA**2   # Hartree units (m c^2 = 1/alpha^2)
s5 = (somm(0, 2, MU) - somm(1, 1, MU)) * HARTREE_HZ / 1e6
V += [('V5 Sommerfeld k = j + 1/2 (mu), bez a_e', s5), ('V5b Sommerfeld + 2 a_e x (LS Pauli, QM)', s5 + split(r3_qm, 2 + 2 * AE) - split(r3_qm, 2.0))]
pr(f"pomiar 2P3/2 - 2P1/2: {units(MEAS)}  (9 911,200(12) Hagley-Pipkin 1994 + 1 057,845(9) Lundeen-Pipkin 1981)")
for name, v in V: pr(f"  {name}: {units(v)};  - pomiar {v-MEAS:+.3f} MHz ({(v/MEAS-1)*100:+.4f} %)")
d = dict(V)
pr(f"H761a (V1 do 1 %): {'TAK' if abs(d['V1 orbita EBK, g = 2']/MEAS-1) < 0.01 else 'NIE'}")
pr(f"H761b (V3: <r^-3> do 1 % QM i rozszczepienie do 1 % pomiaru): {'TAK' if abs(r3_ebk/r3_qm-1) < 0.01 and abs(d['V3 fala EBK 733, g z a_e']/MEAS-1) < 0.01 else 'NIE'}")
pr(f"H761c (V4 do 0,1 %): {'TAK' if abs(d['V4 QM <r^-3>, g z a_e']/MEAS-1) < 1e-3 else 'NIE'}")
# --- Ps (separately) ---
pr("=== Ps (mu = m/2) ===")
pr("  p-Ps 2^1P_1: S = 0 -> brak precesji spin-orbita, jeden poziom (rozszczepienie LS = 0)")
for name, r3 in [('orbita EBK', r3_orbit), ('fala EBK 733', r3_ebk), ('QM', r3_qm)]:
    for gl, g in [('g = 2', 2.0), ('a_e', 2 * (1 + AE))]:
        xi = (3 / 16) * ALPHA**2 * r3 * ((2 * g - 1) / 3) * HARTREE_HZ / 1e6     # MHz per unit L.S ((2g-1)/2 = 3/2 at g = 2)
        pr(f"  o-Ps 2^3P_J, czesc LS, <r^-3> {name}, {gl}: xi = {xi:.3f} MHz; J2 - J1 = 2 xi = {units(2*xi)}; J1 - J0 = xi = {units(xi)}")
nu = {2: 8626.71, 1: 13012.41, 0: 18498.25}; LSq = {2: 1, 1: -1, 0: -2}; T = {2: -0.4, 1: 2.0, 0: -4.0}
A = np.array([[1, LSq[J], T[J]] for J in (2, 1, 0)]); c, xi_th, tau = np.linalg.solve(A, np.array([-nu[J] for J in (2, 1, 0)]))   # 2^3S_1 above 2^3P_J: E_J = E_S - nu_J
pr(f"  teoria QED nu_J (Czarnecki, Melnikov, Yelkhovsky 1999): rozklad E_J = c + xi <L.S> + tau <S12>: xi = {xi_th:.3f} MHz, tau = {tau:.3f} MHz "
   f"(pomiary: nu_0 = 18 501,02(61) MHz, Gurung i in. 2020; nu_2 = 8 627,94(95) MHz, UCL)")
# --- part C: spin precession angle per radial period ---
ang = math.pi * ALPHA**2 / L**2
pr(f"czesc C: kat precesji spinu na okres radialny (orbita EBK n 2, L 1,5, g = 2): pi alfa^2/L^2 = {ang:.4e} rad = {ang/math.pi:.3e} pi "
   f"(poslizg Maslowa -pi: stosunek {ang/math.pi:.2e}); Ps: ten sam rzad (x (3/2) z dwoch spinow): {1.5*ang:.3e} rad")
# --- part B: F basin under classical precession ---
xi = d['V4 QM <r^-3>, g z a_e'] / 1.5; A12, A32 = 59.17, 23.65
aL = (A12 + 2 * A32 / 2 * 1) / 2; aL = (A12 + A32) / 2 if False else None
# solve (4/3) aL - (1/3) aS = A_1/2 ; (2/3) aL + (1/3) aS = A_3/2 = 23.65/2
aL, aS = np.linalg.solve(np.array([[4 / 3, -1 / 3], [2 / 3, 1 / 3]]), np.array([A12, A32 / 2]))
pr(f"czesc B: xi = {xi:.3f} MHz, a_L = {aL:.3f} MHz, a_S = {aS:.3f} MHz (z HFS 2P1/2 59,17 i 2P3/2 23,65 MHz w obrazie wektorowym)")
rng = np.random.default_rng(761)
def rot(v, ax, th):
    ax = ax / np.linalg.norm(ax); return v * np.cos(th) + np.cross(ax, v) * np.sin(th) + ax * np.dot(ax, v) * (1 - np.cos(th))
starts = []
for Jm in (1.0, 2.0):
    for Fm in ((0.5, 1.5) if Jm == 1.0 else (1.5, 2.5)):
        for k in range(20):
            Lm, Sm, Im = 1.5, 1.0, 1.0
            cLS = (Jm * Jm - Lm * Lm - Sm * Sm) / (2 * Lm * Sm); z = np.array([0, 0, 1.0])
            # L and S with resultant J along z
            thL = math.acos((Jm * Jm + Lm * Lm - Sm * Sm) / (2 * Jm * Lm)); Lv = Lm * np.array([math.sin(thL), 0, math.cos(thL)]); Sv = Jm * z - Lv
            ph = rng.uniform(0, 2 * np.pi); Lv, Sv = rot(Lv, z, ph), rot(Sv, z, ph)
            thI = math.acos((Fm * Fm - Jm * Jm - Im * Im) / (2 * Jm * Im)); Iv = Im * np.array([math.sin(thI), 0, math.cos(thI)]); Iv = rot(Iv, z, rng.uniform(0, 2 * np.pi))
            starts.append((Jm, Fm, Lv, Sv, Iv))
Lv = np.array([s[2] for s in starts]); Sv = np.array([s[3] for s in starts]); Iv = np.array([s[4] for s in starts])
J0 = np.array([s[0] for s in starts]); F0 = np.array([s[1] for s in starts])
tw = 2 * np.pi
def f(Lv, Sv, Iv):
    oL = xi * Sv + aL * Iv; oS = xi * Lv + aS * Iv; oI = aL * Lv + aS * Sv
    return tw * np.cross(oL, Lv), tw * np.cross(oS, Sv), tw * np.cross(oI, Iv)
# run 1 had Tso/40 (length drift 1.4e-2), see uzupelnienie
Tso = 1 / (xi * 2.0); dt = Tso / int(os.environ.get('NSUB761', '160')); nper = 10000; steps = int(nper * Tso / dt)
Fn0 = np.linalg.norm(Lv + Sv + Iv, axis=1); Jmax = np.zeros(len(starts)); Fdev = 0.0
for s in range(steps):
    k1 = f(Lv, Sv, Iv); k2 = f(*(v + dt / 2 * k for v, k in zip((Lv, Sv, Iv), k1))); k3 = f(*(v + dt / 2 * k for v, k in zip((Lv, Sv, Iv), k2)))
    k4 = f(*(v + dt * k for v, k in zip((Lv, Sv, Iv), k3)))
    Lv, Sv, Iv = (v + dt / 6 * (a1 + 2 * a2 + 2 * a3 + a4) for v, a1, a2, a3, a4 in zip((Lv, Sv, Iv), k1, k2, k3, k4))
    if s % 50 == 0:
        Jn = np.linalg.norm(Lv + Sv, axis=1); Jmax = np.maximum(Jmax, np.abs(Jn - J0))
        Fdev = max(Fdev, np.abs(np.linalg.norm(Lv + Sv + Iv, axis=1) / Fn0 - 1).max())
jb = (np.linalg.norm(Lv + Sv, axis=1) > 1.5) == (J0 > 1.5)   # basin by nearest EBK length (1 or 2); run 1 used np.round (half to even -- wrong)
pr(f"  {len(starts)} startow, {nper} okresow spin-orbita: max ||J| - |J|_0| = {Jmax.max():.4e} (j 1/2: {Jmax[J0==1].max():.3e}, j 3/2: {Jmax[J0==2].max():.3e}); "
   f"max |F|/|F|_0 - 1 = {Fdev:.2e}; basen j zachowany {jb.sum()}/{len(starts)}; dlugosci |L|, |S|, |I| dryf {np.abs(np.linalg.norm(Lv,axis=1)-1.5).max():.1e}")
pr(f"H761d: {'TAK' if Fdev < 1e-9 and jb.all() and Jmax.max() < 0.05 else 'NIE'}")
open(os.path.join(D, 'out', 'ana761.txt'), 'w').write('\n'.join(out) + '\n')

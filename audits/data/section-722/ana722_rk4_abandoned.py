"""Audit 722 (offline): the clock resonance of 715 (R0: 2<T> = j hbar Omega_r, NR limit of 2D = j hbar Omega_r) on Kepler, the isotropic oscillator and
Coulomb + beta/r^2.  Units mu = k = hbar = omega = 1.  Orbits: 2D RK4 in the Sundman variable s (dt/ds = r), one radial period from periapsis; the
return to periapsis is located by secant steps on v_r from the last full step.  Integrated with the orbit: t, int T dt, int v_r^2 dt (= oint p_r dr).
Step: one radial period in NSTEP steps (step size from the analytic s-period); check with 2 NSTEP."""
import math, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
NSTEP = 4000

def make(kind, beta=0.0):
    if kind == 'osc':
        return (lambda r: 0.5*r*r), (lambda r: r)            # V, V'
    return (lambda r: -1.0/r + beta/(r*r)), (lambda r: 1.0/(r*r) - 2.0*beta/r**3)

def deriv(y, Vp):
    x, z, vx, vz = y[0], y[1], y[2], y[3]; r = math.hypot(x, z); a = -Vp(r)/r
    vr2 = ((x*vx + z*vz)/r)**2
    return [r*vx, r*vz, r*a*x, r*a*z, r, r*0.5*(vx*vx + vz*vz), r*vr2]

def rk4(y, h, Vp):
    k1 = deriv(y, Vp); k2 = deriv([a + 0.5*h*b for a, b in zip(y, k1)], Vp)
    k3 = deriv([a + 0.5*h*b for a, b in zip(y, k2)], Vp); k4 = deriv([a + h*b for a, b in zip(y, k3)], Vp)
    return [a + h/6*(b + 2*c + 2*d + e) for a, b, c, d, e in zip(y, k1, k2, k3, k4)]

def vr(y): return (y[0]*y[2] + y[1]*y[3])/math.hypot(y[0], y[1])

def periapsis(E, L, V):
    lo, hi = 1e-12, 1e3          # smallest root of E - V(r) - L^2/(2r^2) = 0 (bisection on the effective potential's inner wall)
    g = lambda r: E - V(r) - L*L/(2*r*r)
    hi = 1e-6
    while g(hi) < 0: hi *= 1.5
    lo = hi/1.5
    for _ in range(200):
        m = 0.5*(lo + hi)
        if g(m) < 0: lo = m
        else: hi = m
    return 0.5*(lo + hi)

def orbit(E, L, kind, beta=0.0, nstep=NSTEP, s_est=None):
    V, Vp = make(kind, beta); rp = periapsis(E, L, V)
    y = [rp, 0.0, 0.0, L/rp, 0.0, 0.0, 0.0]
    if s_est is None:            # step size only: Coulomb + beta/r^2 is Kepler in r (L -> Lambda), s-period = 2 pi sqrt(a), a = 1/(2|E|);
        if kind == 'osc':         # oscillator: rough pi/sqrt(r_p r_a); the loop below finds the actual return either way
            ra = math.sqrt(E + math.sqrt(max(E*E - L*L, 0.0))); s_est = math.pi/math.sqrt(rp*ra)
        else: s_est = 2*math.pi*math.sqrt(1.0/(2*abs(E)))
    h = s_est/nstep; passed = False; k = 0
    while True:
        yn = rk4(y, h, Vp); v = vr(yn)
        if v < -1e-14: passed = True
        if passed and v >= 0: break
        y = yn; k += 1
        if k > 4*nstep: raise RuntimeError('no return')
    a, b = 0.0, h; fa, fb = vr(y), v          # secant / bisection on the partial step
    for _ in range(60):
        m = b - fb*(b - a)/(fb - fa) if fb != fa else 0.5*(a + b)
        if not (a < m < b): m = 0.5*(a + b)
        fm = vr(rk4(y, m, Vp))
        if abs(fm) < 1e-15: a = b = m; break
        if fm < 0: a, fa = m, fm
        else: b, fb = m, fm
        if b - a < 1e-16*h: break
    yf = rk4(y, 0.5*(a + b), Vp)
    T = yf[4]; return dict(E=E, L=L, Tmean=yf[5]/T, Om=2*math.pi/T, Jr=yf[6]/(2*math.pi), s=k*h + 0.5*(a + b), rp=rp)

out = []; ok = {}
# H722a: Kepler identity and R0 roots
out.append('== Kepler (beta = 0): identity J Omega_r / (2|E|) and 2<T>/(-2E)')
dev = 0.0
for L in (0.5, 1.0, 1.5):
    for e in (0.1, 0.5, 0.8, 0.95):
        N = L/math.sqrt(1 - e*e); E = -0.5/(N*N); o = orbit(E, L, 'kep')
        J = o['Jr'] + L; r1 = J*o['Om']/(2*abs(E)); r2 = 2*o['Tmean']/(-2*E)
        o2 = orbit(E, L, 'kep', nstep=2*NSTEP, s_est=o['s'])
        dev = max(dev, abs(r1 - 1), abs(r2 - 1))
        out.append('  L %.1f e %.2f: J Omega/(2|E|) - 1 %+.2e, <T>/(-E) - 1 %+.2e, J - N %+.2e; 2x steps: J Omega/(2|E|) change %.1e'
                   % (L, e, r1 - 1, r2 - 1, J - N, abs((o2['Jr'] + L)*o2['Om']/(2*abs(E)) - r1)))
def root(L, j, kind, beta):
    Lam = math.sqrt(L*L + 2*beta); lo, hi = Lam*(1 + 1e-7), j + 3.0
    f = lambda N: (lambda o: 2*o['Tmean'] - j*o['Om'])(orbit(-0.5/(N*N), L, kind, beta))
    flo, fhi = f(lo), f(hi)
    for _ in range(45):
        m = 0.5*(lo + hi); fm = f(m)
        if (fm > 0) == (flo > 0): lo, flo = m, fm
        else: hi, fhi = m, fm
    return 0.5*(lo + hi)
rdev = 0.0
for j in (1, 2, 3):
    for L in (0.5, 1.0, 1.5):
        if L >= j: out.append('  R0 j %d L %.1f: no orbit (N = j <= L)' % (j, L)); continue
        N = root(L, j, 'kep', 0.0); o = orbit(-0.5/(N*N), L, 'kep'); J = o['Jr'] + L; rdev = max(rdev, abs(J - j))
        out.append('  R0 j %d L %.1f: root N %.9f, J = J_r + L %.9f (J - j %+.1e), E %.9f' % (j, L, N, J, J - j, -0.5/(N*N)))
ok['H722a'] = dev < 1e-8 and rdev < 1e-6
out.append('H722a (identity to 1e-8: max %.1e; R0 root J = j to 1e-6: max %.1e): %s' % (dev, rdev, ok['H722a']))
# H722b: oscillator
out.append('== isotropic oscillator V = r^2/2 (omega = 1): 2<T>/E and Omega_r')
dmax = 0.0
for E in (1.0, 2.0, 3.7):
    for L in (0.5, 1.0, 1.5):
        if L >= E: continue
        o = orbit(E, L, 'osc'); r1 = 2*o['Tmean']/E; dmax = max(dmax, abs(r1 - 1), abs(o['Om'] - 2))
        out.append('  E %.1f L %.1f: 2<T>/E - 1 %+.2e, Omega_r - 2 %+.2e, J_r %.9f (EBK (E - L)/2 %.9f)' % (E, L, r1 - 1, o['Om'] - 2, o['Jr'], (E - L)/2))
ok['H722b'] = dmax < 1e-8
out.append('  => R0 condition 2<T> = j Omega_r reads E = 2 j hbar omega for every L; with j_phi: E = N hbar omega; QM (N + 3/2) hbar omega:'
           ' the clock misses 3/2 hbar omega (all Maslov phases)')
out.append('H722b (2<T> = E and Omega_r = 2 omega to 1e-8: max %.1e): %s' % (dmax, ok['H722b']))
# H722c: Coulomb + beta/r^2
out.append('== Coulomb + beta/r^2: clock root vs N = j + 2 beta/Lambda; shift vs QM N = j - L + Lambda (n = j, l + 1/2 = L)')
cdev = 0.0; rat01 = []
for beta in (0.01, 0.05):
    for j in (1, 2):
        for L in (0.5, 1.0, 1.5):
            if L >= j: continue
            Lam = math.sqrt(L*L + 2*beta); N = root(L, j, 'kep', beta); pred = j + 2*beta/Lam; cdev = max(cdev, abs(N - pred))
            dq = Lam - L; dc = N - j; rat = dc/dq
            if beta == 0.01: rat01.append(rat)
            out.append('  beta %.2f j %d L %.1f: clock N %.9f, analytic %.9f (diff %+.1e); clock shift %.6f, QM shift %.6f, ratio %.4f;'
                       ' E clock %.7f vs QM %.7f' % (beta, j, L, N, pred, N - pred, dc, dq, rat, -0.5/N**2, -0.5/(j - L + Lam)**2))
ok['H722c'] = cdev < 1e-6 and all(abs(r - 2) < 0.04 for r in rat01)
out.append('H722c (root = analytic to 1e-6: max %.1e; ratio at beta 0.01 within 2 %% of 2: %s): %s'
           % (cdev, ', '.join('%.4f' % r for r in rat01), ok['H722c']))
text = '\n'.join(out); print(text)
open(os.path.join(HERE, 'out', 'ana722.txt'), 'w').write(text + '\n')

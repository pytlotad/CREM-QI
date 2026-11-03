# Audit 691 (offline): a clock of rest frequency (mc^2 + q phi)/hbar on a relativistic Coulomb orbit -- does the phase after an orbit give
# J = n hbar with Maslov and no fitted parameters?  Units m = c = hbar = 1, k = alpha; test electron; L = (l + 1/2), I_r = n_r + 1/2.
# Variants of the wave frequency E' kept in harmony with the clock along the path (E' - p.v = clock rate in LAB, on time average):
#   dt   : the gauge phase S/hbar (q phi dt, audit 603)             -> E' = E
#   dtau : the literal postulate (mc^2 + q phi) dtau / hbar         -> E' = E + <(a/r)(1 - 1/gamma)>
#   m    : the bare clock mc^2 dtau / hbar (audits 380, 603)        -> E' = E + <a/r>
# usage: python3 ana691.py
import math, os
import numpy as np
from scipy.optimize import brentq

A = 7.2973525693e-3                     # CODATA 2018
MEC2 = 510998.95000                     # eV, CODATA 2018
MPME = 1836.15267343                    # CODATA 2018
H_EVS = 4.135667696e-15                 # eV s
C = 299792458.0
MU = MEC2*MPME/(1 + MPME)               # hydrogen reduced mass, eV
U, W = np.polynomial.legendre.leggauss(400)
U = 0.5*math.pi*(U + 1); W = 0.5*math.pi*W   # nodes on [0, pi]


def turning(E, L):
    a2 = E*E - 1; b = 2*E*A; c = A*A - L*L; d = math.sqrt(max(0.0, b*b - 4*a2*c))
    return sorted(((-b + d)/(2*a2), (-b - d)/(2*a2)))


def orbit(E, L):
    r1, r2 = turning(E, L); c0 = 0.5*(r1 + r2); h = 0.5*(r2 - r1)
    r = c0 - h*np.cos(U); s = math.sqrt(1 - E*E)
    ir = float(np.sum(W*s*h*h*np.sin(U)**2/r))/math.pi      # (1/2pi) oint p_r dr
    dt = (E*r + A)/s                                       # dt/du
    return r, dt, ir


def mean(f, E, L):
    r, dt, _ = orbit(E, L); return float(np.sum(W*f(r)*dt)/np.sum(W*dt))


def level(nr, L):
    f = lambda b: orbit(1 - b, L)[2] - (nr + 0.5)
    b0 = A*A/(2*(nr + 0.5 + L)**2); bc = 1 - math.sqrt(1 - A*A/(L*L))     # bc: circular orbit, I_r = 0
    return 1 - brentq(f, 0.5*b0, bc*(1 - 1e-13), xtol=1e-22, rtol=1e-15, maxiter=200)


def kg(nr, l): return (1 + A*A/(nr + 0.5 + math.sqrt((l + 0.5)**2 - A*A))**2)**-0.5


out = ['units: E in mc^2; binding b = 1 - E; a^4 unit mc^2 a^4 = %.6e eV (hydrogen mu c^2 a^4 = %.6e eV)' % (MEC2*A**4, MU*A**4)]
lev = {}; worst = 0.0; worstb = 0.0
for n in (1, 2, 3):
    for l in range(n):
        nr = n - l - 1; L = l + 0.5; E = level(nr, L); Ek = kg(nr, l)
        worst = max(worst, abs(E - Ek)/Ek); worstb = max(worstb, abs((1 - E) - (1 - Ek))/(1 - Ek))
        g = lambda r: E + A/r
        d_tau = mean(lambda r: (A/r)*(1 - 1/g(r)), E, L); d_m = mean(lambda r: A/r, E, L)
        an = A**4/n**4*(n/L - 0.5)
        lev[(n, l)] = dict(dt=E, dtau=E + d_tau, m=E + d_m, an=an, dtau_d=d_tau)
        out.append(f'n {n} l {l}: b_dt {1 - E:.12e} (KG {1 - Ek:.12e})  E_dtau - E_dt {d_tau/A**4:+.9f} a^4 (analytic {an/A**4:+.9f}, rel {d_tau/an - 1:+.2e})'
                   f'  E_m - 1 {E + d_m - 1:+.6e} (= +{(E + d_m - 1)/(1 - E):.6f} x binding)')
out.append(f'H691a: dt = Klein-Gordon, max relative |E - E_KG| {worst:.2e} (binding {worstb:.2e}): {worst <= 1e-12}')
rb = max(abs(v['dtau_d']/v['an'] - 1) for v in lev.values())
out.append(f'H691b: dtau shift = (mc^2 a^4/n^4)(n/(l+1/2) - 1/2) to relative {rb:.2e} (order a^2 = {A*A:.1e}): {rb <= 1e-3}')


def line(lab, dE):   # dE in mc^2 -> hydrogen frequency, wavelength, eV, ueV
    e = dE*MU; f = e/H_EVS
    return f'{lab} {f/1e6:+12.3f} MHz  lambda {C/abs(f)*1e2:9.4f} cm  {e:+.6e} eV  {e*1e6:+10.4f} ueV'


out.append('hydrogen intervals (mu c^2 scaling, recoil a^4 m/M neglected); Dirac centroid (j-weighted) = KG at a^4 for l >= 1')
for v in ('dt', 'dtau', 'm'):
    out.append(line(f'  {v:4s} 2P - 2S         ', lev[(2, 1)][v] - lev[(2, 0)][v]))
    out.append(line(f'  {v:4s} 3D - 3P         ', lev[(3, 2)][v] - lev[(3, 1)][v]))
dir_ = lambda n, j: (1 + A*A/(n - j - 0.5 + math.sqrt((j + 0.5)**2 - A*A))**2)**-0.5
c2p = (2*dir_(2, 0.5) + 4*dir_(2, 1.5))/6; c3p = (2*dir_(3, 0.5) + 4*dir_(3, 1.5))/6; c3d = (4*dir_(3, 1.5) + 6*dir_(3, 2.5))/10
out.append(line('  Dirac 2P_c - 2S1/2  ', c2p - dir_(2, 0.5)))
out.append(line('  Dirac 3D_c - 3P_c   ', c3d - c3p))
meas = (2*(-1057.845) + 4*9911.200)/6
out.append(f'  measured 2P_c - 2S1/2 {meas:+.2f} MHz (Lundeen & Pipkin 1981 2S1/2 - 2P1/2 1057.845 MHz; Hagley & Pipkin 1994 2P3/2 - 2S1/2 9911.200 MHz)'
           f'  lambda {C/(meas*1e6)*1e2:.4f} cm  {meas*1e6*H_EVS:.6e} eV  {meas*1e6*H_EVS*1e6:.4f} ueV')
s32 = (lev[(3, 2)]['dtau'] - lev[(3, 1)]['dtau'])*(lev[(3, 2)]['dt'] - lev[(3, 1)]['dt'])
s21 = (lev[(2, 1)]['dtau'] - lev[(2, 0)]['dtau'])*(lev[(2, 1)]['dt'] - lev[(2, 0)]['dt'])
out.append(f'H691c: dtau reverses the sign of the l ordering against KG / Dirac centroid (3D - 3P, 2P - 2S): {s32 < 0 and s21 < 0}')
mm = all(v['m'] > 1 for v in lev.values())
out.append(f'H691d: bare clock m: E\' > mc^2 for every level (levels upside down, no Bohr condition): {mm}')

# H691e: retarded wave at speed c, nonrelativistic Kepler orbit n = 2, e = 0.5 (atomic units, c = 1/alpha)
cA = 1/A; a = 4.0; e = 0.5; T = 2*math.pi*a**1.5


def pos(t):
    M = 2*math.pi*t/T; Ea = M
    for _ in range(60): Ea -= (Ea - e*math.sin(Ea) - M)/(1 - e*math.cos(Ea))
    return np.array([a*(math.cos(Ea) - e), a*math.sqrt(1 - e*e)*math.sin(Ea)])


t0 = 10*T; one = 0; delays = []
tr = np.linspace(t0 - 10*T, t0, 200001); P = np.array([pos(x) for x in tr])
for j in range(100):
    xf = pos(t0 + T*(j + 0.5)/100)                                   # a point of the orbit away from the particle
    g = cA*(t0 - tr) - np.linalg.norm(P - xf, axis=1)
    roots = np.sum(np.sign(g[1:]) != np.sign(g[:-1])); one += roots == 1
    k = int(np.argmin(np.abs(g))); delays.append((t0 - tr[k])/T)
out.append(f'H691e: retarded time unique over 10 previous orbits in {one}/100 field points; delay {min(delays):.2e}-{max(delays):.2e} of a period: {one == 100}')

# H691f: the clock rate (dtau variant) on both sides of the turning points (function of r only, so continuous): no -pi/2 jump
E = lev[(2, 1)]['dt']; r1, r2 = turning(E, 1.5); rate = lambda r: (1 - A/r)/(E + A/r)
jump = max(abs(rate(r*(1 + 1e-9)) - rate(r*(1 - 1e-9))) for r in (r1, r2))
out.append(f'H691f: clock rate jump across the turning points {jump:.1e} (smooth; Maslov -pi/2 not from the clock): {jump < 1e-9}')
text = '\n'.join(out); print(text)
os.makedirs('out', exist_ok=True); open('out/ana691.txt', 'w').write(text + '\n')

"""Audit 722 (offline): the identity J Omega = 2|E| (Kepler) and the clock resonance R0 (2 <T> = j hbar Omega_r, the non-relativistic limit of
2 D = j hbar Omega_r from 715) for non-Keplerian potentials.  Units mu = k = hbar = omega = 1.
Radial motion only (the R0 condition does not involve the angle): r'' = L^2/r^3 - V'(r), from periapsis (p_r = 0) to apoapsis (event p_r = 0)
= half a radial period, by symmetry.  Quadratures carried as ODE states: t, int T dt, int p_r^2 dt.
  <T> = int T dt / (T_r/2),  Omega_r = 2 pi / T_r,  J_r = (1/2 pi) oint p_r dr = (1/pi) int_half p_r^2 dt.
Integrator: scipy DOP853, rtol 1e-13, atol 1e-15; control: the same with rtol 1e-11 (difference printed).
(The first version -- pure-Python RK4, as written in prereg722 -- produced no result in ~3 h and was stopped; it is kept as ana722_rk4_abandoned.py.)
"""
import math, os, sys
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__))


def pot(kind, beta):
    if kind == 'kep': return (lambda r: -1/r), (lambda r: 1/r**2)
    if kind == 'cb': return (lambda r: -1/r + beta/r**2), (lambda r: 1/r**2 - 2*beta/r**3)
    if kind == 'osc': return (lambda r: 0.5*r*r), (lambda r: r)


def periapsis(kind, E, L, beta=0.0):
    if kind == 'osc': return math.sqrt(E - math.sqrt(E*E - L*L))
    Lam2 = L*L + 2*beta                     # E r^2 + r - Lam^2/2 = 0, smaller root
    return (1 - math.sqrt(1 + 2*E*Lam2))/(-2*E)


def orbit(kind, E, L, beta=0.0, rtol=1e-13):
    V, dV = pot(kind, beta); r0 = periapsis(kind, E, L, beta)
    def f(t, y):
        r, pr = y[0], y[1]
        return [pr, L*L/r**3 - dV(r), 1.0, 0.5*(pr*pr + L*L/(r*r)), pr*pr]
    def apo(t, y): return y[1]
    apo.terminal = True; apo.direction = -1
    # scale of the half period for t_max
    tmax = 50*math.pi*(1/(2*abs(E)))**1.5 if kind != 'osc' else 10*math.pi
    # start a hair after periapsis so the event p_r = 0 at t = 0 is not caught
    sol = solve_ivp(f, (0, tmax), [r0, 0.0, 0.0, 0.0, 0.0], method='DOP853', rtol=rtol, atol=1e-15, events=apo,
                    first_step=None, dense_output=False)
    y = sol.y_events[0][-1] if len(sol.t_events[0]) else None
    if y is None: raise RuntimeError('no apoapsis %s E %g L %g' % (kind, E, L))
    half = y[2]; Tr = 2*half
    return {'E': E, 'L': L, 'T': y[3]/half, 'Om': 2*math.pi/Tr, 'Jr': y[4]/math.pi, 'ra': y[0], 'rp': r0,
            'Eerr': 0.5*y[1]**2 + L*L/(2*y[0]**2) + V(y[0]) - E}


out = []
def P(s): out.append(s); print(s, flush=True)

P('Audit 722: J Omega = 2|E| and the clock condition 2<T> = j Omega_r (mu = k = hbar = omega = 1); DOP853 rtol 1e-13, control rtol 1e-11')
# ---- H722a: Kepler identity on a grid, and R0 roots
P('\n(A) Kepler: identity J Omega_r / (2|E|) on a grid (J = J_r + L)')
worst = 0.0; worstc = 0.0
for E in (-0.5, -0.125, -0.05):
    Lmax = 1/math.sqrt(-2*E)
    for e in (0.1, 0.5, 0.8, 0.95):
        L = Lmax*math.sqrt(1 - e*e)
        o = orbit('kep', E, L); oc = orbit('kep', E, L, rtol=1e-11)
        J = o['Jr'] + L; idv = J*o['Om']/(2*abs(E)); vir = o['T']/(-E)
        worst = max(worst, abs(idv - 1), abs(vir - 1)); worstc = max(worstc, abs((oc['Jr'] + L)*oc['Om']/(2*abs(E)) - idv))
        P('  E %+.4f e %.2f L %.6f: J %.10f  J Om/2|E| - 1 %+.2e  <T>/|E| - 1 %+.2e  (energy drift %.1e)' % (E, e, L, J, idv - 1, vir - 1, o['Eerr']))
P('  max |identity - 1|, |virial - 1| = %.2e; control (rtol 1e-11) changes identity by <= %.1e' % (worst, worstc))
okA1 = worst < 1e-8
P('  R0 roots (2<T> = j Omega_r at fixed L), Kepler: predicted J = j exactly, E = -1/(2 j^2)')
worstr = 0.0
for j in (1, 2, 3):
    for L in (0.5, 1.0, 1.5):
        if L >= j: P('  j %d L %.1f: skipped (no Kepler orbit with J = j and L >= j)' % (j, L)); continue
        g = lambda E: 2*orbit('kep', E, L)['T'] - j*orbit('kep', E, L)['Om']
        Elo = -1/(2*L*L)*(1 - 1e-6); Ehi = -1e-4
        E0 = brentq(g, Elo, Ehi, xtol=1e-15, rtol=1e-14)
        o = orbit('kep', E0, L); J = o['Jr'] + L; worstr = max(worstr, abs(J - j))
        P('  j %d L %.1f: root E %.12f (Bohr %.12f), J = %.10f, J - j %+.2e' % (j, L, E0, -1/(2*j*j), J, J - j))
okA2 = worstr < 1e-6
P('H722a (identity to 1e-8 and R0 root at J = j to 1e-6): %s (identity %.1e, root %.1e)' % (okA1 and okA2, worst, worstr))

# ---- H722b: oscillator
P('\n(B) Isotropic oscillator V = r^2/2: 2<T> = E and Omega_r = 2 for every (E, L)')
worstb = 0.0
for E in (1.0, 2.5, 4.0, 7.0):
    for frac in (0.05, 0.3, 0.6, 0.9):
        L = frac*E; o = orbit('osc', E, L)
        d1 = 2*o['T']/E - 1; d2 = o['Om']/2 - 1; worstb = max(worstb, abs(d1), abs(d2))
        Jr = o['Jr']
        P('  E %.2f L %.3f: 2<T>/E - 1 %+.2e  Omega_r/2 - 1 %+.2e  J_r %.10f (EBK: E = 2 J_r + L -> %.10f)' % (E, L, d1, d2, Jr, 2*Jr + L))
okB = worstb < 1e-8
P('  => clock R0 condition 2<T> = j Omega_r <=> E = 2 j: E_clock = N hbar omega (N = 2j, even), any L')
P('     QM (= EBK with Maslov): E = (N + 3/2) hbar omega -- clock spectrum lacks the zero point 3/2 hbar omega (N = 2: 2 vs 3.5)')
P('H722b (2<T> = E and Omega_r = 2 to 1e-8): %s (max %.1e)' % (okB, worstb))

# ---- H722c: Coulomb + beta/r^2
P('\n(C) Coulomb + beta/r^2: clock root N vs prediction j + 2 beta/Lambda; shift ratio clock/QM')
worstc_ = 0.0; ratios = {}
for beta in (0.01, 0.05):
    for j in (1, 2):
        for L in (0.5, 1.0, 1.5):
            if L >= j: continue
            Lam = math.sqrt(L*L + 2*beta)
            g = lambda E: 2*orbit('cb', E, L, beta)['T'] - j*orbit('cb', E, L, beta)['Om']
            Emin = -1/(2*Lam*Lam)*(1 - 1e-6)
            E0 = brentq(g, Emin, -1e-4, xtol=1e-15, rtol=1e-14)
            N = 1/math.sqrt(-2*E0); Npred = j + 2*beta/Lam; worstc_ = max(worstc_, abs(N - Npred))
            o = orbit('cb', E0, L, beta); Ncheck = o['Jr'] + Lam
            NQM = (j - L) + Lam                  # n_r + 1/2 = j - L, sqrt((l+1/2)^2 + 2 beta) = Lambda
            sc = N - j; sq = NQM - j; r = sc/sq; ratios[(beta, j, L)] = r
            P('  beta %.2f j %d L %.1f: clock N %.10f (pred %.10f, diff %+.1e; J_r + Lambda %.10f)  QM N %.10f  shift clock %+.6f  QM %+.6f  ratio %.5f'
              % (beta, j, L, N, Npred, N - Npred, Ncheck, NQM, sc, sq, r))
            # energies, a.u. (Ha)
            P('      E clock %.10f  E QM %.10f  E Kepler(j) %.10f Ha' % (-1/(2*N*N), -1/(2*NQM*NQM), -1/(2*j*j)))
okC1 = worstc_ < 1e-6
okC2 = all(abs(r - 2) < 0.04 for (b, j, L), r in ratios.items() if b == 0.01)
P('H722c (root = j + 2 beta/Lambda to 1e-6; ratio within 2 %% of 2 at beta = 0.01): %s (root %.1e; ratios beta 0.01: %s)'
  % (okC1 and okC2, worstc_, ', '.join('%.4f' % r for (b, j, L), r in ratios.items() if b == 0.01)))
os.makedirs(os.path.join(HERE, 'out'), exist_ok=True)
open(os.path.join(HERE, 'out', 'ana722.txt'), 'w').write('\n'.join(out) + '\n')

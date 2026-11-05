# Audit 712 (offline): the mechanical Compton clock with a primary-frequency component of the particle's field.
# Part 1: on relativistic Coulomb orbits (H = sum sqrt(m_i^2 c^4 + p^2 c^2) - 1/r, relative coordinates, atomic units) the summed mean clock
#   offsets sum_i <m_i c^2/gamma_i - m_i c^2> equal the binding energy E (virial for 1/r); per-particle share against (m_j/M) E.
# Part 2: hydrogen (infinite nucleus), Kepler orbits; the oscillating component driven by the clock, e^{-i theta}, theta' = -T + s_phi V; with the
#   component coupled to the potential (s_g = 1) its bound modes phi_nlm grow at resonances (E_n - <theta'>) T_p = 2 pi j; peak strength
#   G = sum_lm |<phi*_nlm e^{-i theta_p} e^{i Delta t}>_period|^2.  Variants W0 (theta = 0), W3 (s_phi = 0), W4 (s_phi = 1).
# usage: python3 ana712.py
import os, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import genlaguerre, sph_harm_y, zeta

ALPHA = 7.2973525693e-3; C = 1/ALPHA; MP = 1836.15267343; HARTREE_EV = 27.211386245988
out = []


def part1():
    ok = True
    out.append('PART 1: relativistic Coulomb orbits, mean clock offsets over 20 radial periods (atomic units, c = 1/alpha)')
    for name, m1, m2 in (('hydrogen, infinite nucleus', 1.0, None), ('hydrogen, proton', 1.0, MP), ('positronium', 1.0, 1.0)):
        mu = m1 if m2 is None else m1*m2/(m1 + m2)
        for n in (1, 2, 3):
            for e in (0.001, 0.5, 0.9):
                a = n*n/mu; r0 = a*(1 - e); p0 = math.sqrt(mu*(1 + e)/(a*(1 - e)))   # relative momentum at periapsis (vis-viva, k = 1)
                def en(p2):
                    E1 = math.sqrt(m1*m1*C**4 + p2*C*C); E2 = 0.0 if m2 is None else math.sqrt(m2*m2*C**4 + p2*C*C)
                    return E1, E2
                def rhs(t, y):
                    x, yy, px, py = y[:4]; r = math.hypot(x, yy); p2 = px*px + py*py; E1, E2 = en(p2)
                    k = C*C*(1/E1 + (0.0 if m2 is None else 1/E2))
                    s1 = -p2*C*C*m1*C*C/(E1*(E1 + m1*C*C)); s2 = 0.0 if m2 is None else -p2*C*C*m2*C*C/(E2*(E2 + m2*C*C))   # m c^2/gamma - m c^2 without cancellation
                    return [k*px, k*py, -x/r**3, -yy/r**3, s1, s2]
                def peri(t, y): return y[0]*y[2] + y[1]*y[3]
                peri.direction = 1
                Tp = 2*math.pi*math.sqrt(a**3/mu)
                sol = solve_ivp(rhs, (0, 20.5*Tp), [r0, 0.0, 0.0, p0, 0.0, 0.0], method='DOP853', rtol=1e-13, atol=[1e-14*a, 1e-14*a, 1e-14*p0, 1e-14*p0, 1e-14, 1e-14], events=peri, dense_output=True)
                ev = sol.t_events[0]
                if len(ev) < 3:   # no periapsis events: average over the run
                    ta, tb = 0.0, 20.0*Tp
                else:
                    ta, tb = ev[0], ev[-1]
                ya, yb = sol.sol(ta), sol.sol(tb)
                s1 = (yb[4] - ya[4])/(tb - ta); s2 = (yb[5] - ya[5])/(tb - ta)
                E1, E2 = en(p0*p0); H = p0*p0*C*C/(E1 + m1*C*C) + (0.0 if m2 is None else p0*p0*C*C/(E2 + m2*C*C)) - 1/r0
                rel = abs(s1 + s2 - H)/abs(H); ok &= rel < 1e-9
                share = '' if m2 is None else f' | share particle 1 {s1/H:.9f} vs m2/M {m2/(m1 + m2):.9f} (diff {s1/H - m2/(m1 + m2):+.2e}), particle 2 {s2/H:.9f}'
                out.append(f'  {name:27s} n {n} e {e:5.3f}: E {H:+.12e} Ha, sum of clock offsets {s1 + s2:+.12e}, |sum - E|/|E| {rel:.1e}{share}')
    out.append(f'H712a (|sum of clock offsets - E|/|E| < 1e-9 everywhere): {ok}')


def R(n, l, r):
    norm = math.sqrt((2/n)**3*math.factorial(n - l - 1)/(2*n*math.factorial(n + l)))
    x = 2*r/n
    return norm*np.exp(-r/n)*x**l*genlaguerre(n - l - 1, 2*l + 1)(x)


def strength(nu, L, n, Delta, sphi, M=4000):
    """time average over one period of phi*_nlm e^{-i theta_p} e^{i Delta t}, summed |.|^2 over l, m (orbit in the xy plane)"""
    a = nu*nu; e = math.sqrt(max(0.0, 1 - (L/nu)**2))
    u = (np.arange(M) + 0.5)*2*math.pi/M; r = a*(1 - e*np.cos(u))
    phi = 2*np.arctan2(math.sqrt(1 + e)*np.sin(u/2), math.sqrt(1 - e)*np.cos(u/2))
    t = nu**3*(u - e*np.sin(u)); w = (1 - e*np.cos(u))/M
    thp = {None: 0.0*u, 0: -nu*e*np.sin(u), 1: -2*nu*e*np.sin(u)}[sphi]
    ph = np.exp(-1j*thp + 1j*Delta*t)*w
    G = 0.0
    for l in range(n):
        Rl = R(n, l, r)
        for m in range(-l, l + 1):
            Y = sph_harm_y(l, m, math.pi/2, phi)
            G += abs(np.sum(np.conj(Rl*Y)*ph))**2
    return G


def resonances(W, n, L):
    """(nu, j) with (E_n - <theta'>) T_p = 2 pi j, nu in [max(0.9, L + 0.01), 4.6]"""
    lo, hi = max(0.9, L + 0.01), 4.6
    if W == 'W0': f = lambda nu: -nu**3/(2*n*n)                 # Delta = E_n, T_p = 2 pi nu^3 -> Delta T_p/(2 pi)
    elif W == 'W3': f = lambda nu: 0.5*nu*(1 - nu*nu/(n*n))
    else: f = lambda nu: 0.5*nu*(3 - nu*nu/(n*n))
    grid = np.linspace(lo, hi, 20001); fv = f(grid); res = []
    for j in range(int(math.floor(fv.min())) - 1, int(math.ceil(fv.max())) + 2):
        g = fv - j
        for i in np.nonzero(np.sign(g[:-1])*np.sign(g[1:]) <= 0)[0]:
            if g[i] == 0 and i > 0 and g[i - 1] == 0: continue
            nu = brentq(lambda x: f(x) - j, grid[i], grid[i + 1]) if g[i] != g[i + 1] else grid[i]
            if all(abs(nu - q) > 1e-9 for q, _ in res): res.append((nu, j))
    return sorted(res)


def part2():
    out.append('PART 2: hydrogen, Kepler orbits, resonance strength G of the bound modes of the potential-coupled component (atomic units)')
    okB = okC = okD = True
    for W, sphi, mean in (('W3', 0, lambda nu: -0.5/nu**2), ('W4', 1, lambda nu: -1.5/nu**2), ('W0', None, lambda nu: 0.0)):
        for L in (0.5, 1.0, 1.5):
            for n in (1, 2, 3, 4):
                res = resonances(W, n, L)
                if not res: out.append(f'  {W} L {L} mode n {n}: no resonance in range'); continue
                G = [(strength(nu, L, n, -0.5/n**2 - mean(nu), sphi), nu, j) for nu, j in res]
                G.sort(reverse=True); top = G[0]; j0 = [g for g in G if g[2] == 0]
                second = G[1] if len(G) > 1 else (0.0, float('nan'), 0)
                out.append(f'  {W} L {L} mode n {n}: {len(G)} resonances; largest G {top[0]:.3e} at nu {top[1]:.4f} (j {top[2]}), next {second[0]:.3e} at nu {second[1]:.4f}'
                           f' (j {second[2]})' + (f' | j = 0 at nu {j0[0][1]:.4f}: G {j0[0][0]:.3e}' if j0 else ' | no j = 0 resonance'))
                if W == 'W3' and n > L: okB &= abs(top[1] - n) < 0.01
                if W == 'W4' and n in (1, 2): okC &= abs(top[1] - math.sqrt(3)*n) < 0.01
                if W == 'W0': okD &= not j0 and all(abs(nu - n) > 0.01 for _, nu, jj in G if jj == 0)
    out.append(f'H712b (W3: largest peak of every mode n = 1-4, every L, at |nu - n| < 0.01): {okB}')
    out.append(f'H712c (W4: largest peak of modes n = 1, 2 at nu = sqrt(3) n +- 0.01): {okC}')
    out.append(f'H712d (W0: no j = 0 resonance; peaks only on sidebands nu = (2 n^2 |j|)^(1/3)): {okD}')


def part3():
    z = zeta(1.5); alt = (1 - 2**-0.5)*z
    out.append('PART 3 (analytic): free component (W1, W2): return after k orbits ~ k^(-3/2); sum k^(-3/2) e^{i k Phi}: peak zeta(3/2) = '
               f'{z:.4f}, trough (1 - 2^(-1/2)) zeta(3/2) = {alt:.4f}, contrast {z/alt:.3f}; Phi per orbit W1 = -J/(2 hbar) (peak J = 2 n h), '
               'W2 = -3J/(2 hbar) (peak J = 2/3 n h).')
    ry_mu = 0.25*HARTREE_EV
    out.append(f'  Ps: each clock (m_j/M) E = E/2 (part 1); single-particle field (mass m) in the partner potential: levels -{2*ry_mu:.4f}/n^2 eV ->'
               f' E/2 = level gives E = -{4*ry_mu:.4f}/n^2 eV vs Ps -{ry_mu:.4f}/n^2 eV (factor 4); relative-coordinate field (mu) with the sum of both'
               ' clocks (= E, part 1) gives -Ry(mu)/n^2.')
    e12 = 0.75*13.598287264   # reduced-mass Bohr 1S-2S, eV (Ry_inf m_p/(m_p + m_e))
    out.append(f'  hydrogen 1S-2S at the W3 resonances (Bohr, reduced mass): {e12:.5f} eV = {e12/4.135667696e-15/1e12:.3f} THz = {1239.84198/e12:.4f} nm'
               ' vs measured 2466.061413 THz (Parthey et al., PRL 107, 203001 (2011)); the rest is fine structure and Lamb shift.')


if __name__ == '__main__':
    part1(); part2(); part3()
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana712.txt', 'w').write(text + '\n')

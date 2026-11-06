"""Audit 732 (offline): E1 rates from the matrix element of the EBK wave (Langer + uniform Airy at both turning points, EBK energies) vs exact QM,
the engine and measurement.  Atomic units, infinite nucleus; reduced mass applied to the rates.  See prereg732.txt."""
import math, os
import numpy as np
from scipy.special import airy, genlaguerre, factorial
HERE = os.path.dirname(os.path.abspath(__file__))
out = []
def P(s): out.append(s); print(s, flush=True)
ALPHA = 1/137.035999177; TAU_AU = 2.4188843265864e-17; MU_H = 1/(1 + 1/1836.152673426)
r = np.concatenate([np.geomspace(1e-9, 1e-3, 4000)[:-1], np.linspace(1e-3, 200.0, 2000001)])


def cumtrapz(y, x):
    return np.concatenate([[0.0], np.cumsum(0.5*(y[1:] + y[:-1])*np.diff(x))])


def u_qm(n, l):
    rho = 2*r/n
    R = math.sqrt((2/n)**3*factorial(n - l - 1)/(2*n*factorial(n + l)))*np.exp(-rho/2)*rho**l*genlaguerre(n - l - 1, 2*l + 1)(rho)
    u = r*R
    return u/math.sqrt(np.trapezoid(u*u, r))


def u_ebk(n, l):
    E = -0.5/n**2; Lam2 = (l + 0.5)**2; nr = n - l - 1
    Q = 2*(E + 1/r) - Lam2/r**2
    disc = math.sqrt(1 + 2*E*Lam2); r1 = (-1 + disc)/(2*E) if False else (1 - disc)/(-2*E); r2 = (1 + disc)/(-2*E)
    sq = np.sqrt(np.abs(Q))
    # phase from r1 (inner) and from r2 (outer); signed so that zeta > 0 in the allowed region
    C = cumtrapz(sq, r)
    C1 = np.interp(r1, r, C); C2 = np.interp(r2, r, C)
    ph_in = np.where(r >= r1, C - C1, -(C1 - C))          # >0 allowed (r > r1), <0 forbidden
    ph_out = np.where(r <= r2, C2 - C, -(C - C2))          # >0 allowed (r < r2), <0 forbidden
    def branch(ph):
        zeta = np.sign(ph)*(1.5*np.abs(ph))**(2/3)
        with np.errstate(divide='ignore', invalid='ignore'):
            pref = np.where(np.abs(Q) > 1e-300, (zeta/Q)**0.25, 0.0)
        pref = np.nan_to_num(pref, nan=0.0)
        return pref*airy(-zeta)[0]
    ui, uo = branch(ph_in), branch(ph_out)
    rm = min(max(Lam2, r1 + 0.05*(r2 - r1)), r2 - 0.05*(r2 - r1))
    k = np.searchsorted(r, rm)
    u = np.where(r < rm, ui, (-1)**nr*uo)
    jump = abs(ui[k] - (-1)**nr*uo[k])/max(abs(ui[k]), 1e-30)
    # fix the overall sign to match QM convention at large r (positive tail times (-1)^nr handled by normalization sign below)
    u = u/math.sqrt(np.trapezoid(u*u, r))
    return u, jump, (r1, r2, rm)


def rad(ua, ub): return np.trapezoid(ua*ub*r, r)


def rate(n_i, l_i, n_f, l_f, R, mu=1.0):
    om = 0.5/n_f**2 - 0.5/n_i**2; lg = max(l_i, l_f)
    A = (4/3)*ALPHA**3*om**3*(lg/(2*l_i + 1))*R*R      # a.u., infinite nucleus
    return A*mu/TAU_AU                                   # s^-1, reduced mass (A ~ mu)


P('Audit 732: E1 matrix elements of the EBK wave (Langer + uniform Airy at both turning points, E = -1/2n^2) vs exact QM; grid to r = 200, 2e6 points')
U = {}
for (n, l) in ((1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2), (4, 1), (4, 3)):
    uq = u_qm(n, l); ue, jump, tp = u_ebk(n, l)
    if np.trapezoid(uq*ue, r) < 0: ue = -ue
    U[(n, l)] = (uq, ue)
    P('  state n %d l %d: turning points %.4f / %.4f, switch at %.4f, branch mismatch at switch %.2e, overlap <QM|EBK> = %.6f' %
      (n, l, tp[0], tp[1], tp[2], jump, np.trapezoid(uq*ue, r)))
P('')
trans = [((2, 1), (1, 0)), ((3, 1), (1, 0)), ((3, 1), (2, 0)), ((3, 0), (2, 1)), ((3, 2), (2, 1)), ((4, 1), (1, 0)), ((4, 3), (3, 2))]
okC = True; res = {}
for (i, f) in trans:
    Rq = rad(U[i][0], U[f][0]); Re = rad(U[i][1], U[f][1])
    if Re*Rq < 0: Re = -Re
    res[(i, f)] = (Rq, Re); okC &= abs(Re/Rq - 1) < 0.03
    P('  %d%s -> %d%s: R_QM %+.6f  R_EBK %+.6f  ratio %.5f  (rate ratio %.5f)' % (i[0], 'spdf'[i[1]], f[0], 'spdf'[f[1]], Rq, Re, Re/Rq, (Re/Rq)**2))
Rq21, Re21 = res[((2, 1), (1, 0))]
okA = abs(abs(Rq21) - 2**7*math.sqrt(6)/3**5) < 1e-7
tq = 1/rate(2, 1, 1, 0, Rq21, MU_H); te = 1/rate(2, 1, 1, 0, Re21, MU_H)
okA &= abs(tq*1e9 - 1.5961) < 1e-3
P('\nH732a: R_QM(2p->1s) %.8f vs 2^7 sqrt6/3^5 = %.8f; tau_QM(2p) = %.5f ns (reduced mass): %s' % (abs(Rq21), 2**7*math.sqrt(6)/3**5, tq*1e9, okA))
P('hydrogen 2p lifetime: EBK wave %.5f ns; QM %.5f ns; engine 1.6747 ns; measured 1.600(4) ns (Bickel, Goodman 1966)' % (te*1e9, tq*1e9))
P('  EBK - measured %+.4f ns (%+.1f sigma); engine - measured +0.0747 ns (+18.7 sigma)' % (te*1e9 - 1.600, (te*1e9 - 1.600)/0.004))
okB = abs(te*1e9 - 1.600) < 0.0747
P('H732b (EBK closer to measurement than the engine): %s' % okB)
# 3p
A3q = rate(3, 1, 1, 0, res[((3, 1), (1, 0))][0], MU_H) + rate(3, 1, 2, 0, res[((3, 1), (2, 0))][0], MU_H)
A3e = rate(3, 1, 1, 0, res[((3, 1), (1, 0))][1], MU_H) + rate(3, 1, 2, 0, res[((3, 1), (2, 0))][1], MU_H)
P('hydrogen 3p lifetime: EBK wave %.4f ns; QM %.4f ns; engine 5.5848 ns; measured 5.58(13) ns' % (1e9/A3e, 1e9/A3q))
P('H732c (|R_EBK/R_QM - 1| < 3 %% for all seven transitions): %s' % okC)
# Ps
tPs_e = 2/rate(2, 1, 1, 0, Re21, 1.0); tPs_q = 2/rate(2, 1, 1, 0, Rq21, 1.0)
P('H732d (report) Ps 2P -> 1S (p-Ps 2^1P1 -> 1^1S0 and o-Ps 2^3P_J -> 1^3S1, same radial element in a_Ps, omega/2): EBK %.4f ns, QM %.4f ns' % (tPs_e*1e9, tPs_q*1e9))
lines = []
for (i, f) in trans:
    om = 0.5/f[0]**2 - 0.5/i[0]**2
    lines.append('%d%s->%d%s: %.4f eV = %.3f nm' % (i[0], 'spdf'[i[1]], f[0], 'spdf'[f[1]], om*27.211386245988*MU_H, 1239.8419843320026/(om*27.211386245988*MU_H)))
P('lines (Bohr, reduced mass): ' + '; '.join(lines))
open(os.path.join(HERE, 'out', 'ana732.txt'), 'w').write('\n'.join(out) + '\n')

"""Audit 730 (offline): Maslov indices on the two-spin torus in the Schwinger polarization (spin = two oscillators, configuration space R^2 per spin).
z in C^4 = (a1, a2, b1, b2), x = sqrt2 Re z, p = sqrt2 Im z; quadratic generators H = z^+ M z, flow z(t) = exp(-i M t) z0.  See prereg730.txt."""
import math, os
import numpy as np
from scipy.linalg import expm
HERE = os.path.dirname(os.path.abspath(__file__))
out = []
def P(s): out.append(s); print(s, flush=True)
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex); I2 = np.eye(2); Z2 = np.zeros((2, 2))
def blk(A, B): return np.block([[A, Z2], [Z2, B]])
M_Na = blk(I2, Z2); M_Nb = blk(Z2, I2); M_Sz = 0.5*blk(sz, sz)
def M_Sn(n): s = n[0]*sx + n[1]*sy + n[2]*sz; return 0.5*blk(s, s)


def spinor(direction, length, chi):
    th = math.acos(max(-1, min(1, direction[2]))); ph = math.atan2(direction[1], direction[0])
    return math.sqrt(2*length)*np.exp(1j*chi)*np.array([math.cos(th/2), np.exp(1j*ph)*math.sin(th/2)])


def spins(z):
    a, b = z[:2], z[2:]
    s1 = 0.5*np.real([np.vdot(a, sx @ a), np.vdot(a, sy @ a), np.vdot(a, sz @ a)])
    s2 = 0.5*np.real([np.vdot(b, sx @ b), np.vdot(b, sy @ b), np.vdot(b, sz @ b)])
    return s1, s2


def tangent(z, n):
    cols = [-1j*(M @ z) for M in (M_Na, M_Nb, M_Sn(n), M_Sz)]
    R = np.stack([np.concatenate([math.sqrt(2)*v.real, math.sqrt(2)*v.imag]) for v in cols], 1)   # 8 x 4 real
    Q, _ = np.linalg.qr(R)
    return Q[:4] + 1j*Q[4:]                                                                       # X + i P, 4 x 4 unitary


def cycle(z0, Mg, n, steps=12000):
    T = 2*np.pi; dt = T/steps; Ustep = expm(-1j*Mg*dt)
    z = z0.copy(); ph = []; act = 0.0
    for k in range(steps + 1):
        s1_, s2_ = spins(z); St = s1_ + s2_; nt = St/np.linalg.norm(St)          # Addendum 1: |S| field at the current point
        d = np.linalg.det(tangent(z, nt)); ph.append(np.angle(d*d))
        if k < steps:
            zn = Ustep @ z
            zm = 0.5*(z + zn); zd = (zn - z)/dt
            act += -np.imag(np.vdot(zm, zd))*dt
            z = zn
    mu = -(np.unwrap(np.array(ph))[-1] - ph[0])/(2*np.pi)
    return mu, act, np.abs(z - z0).max()


P('Audit 730 (run 2, Addendum 1: n(t) = S^(z(t))): Maslov indices on the two-spin torus, Schwinger polarization; cycles C_a (N_a), Z (S_z + (N_a + N_b)/2), G (S_n + (N_a + N_b)/2), 2 pi each')
okA = True; muG = []; muZ = []
for s_len in (1.0, 1.5, 10.5):
    for Sfrac in (0.3, 0.6, 0.9):
        Sl = Sfrac*2*s_len; ca = Sl/(2*s_len); sa = math.sqrt(1 - ca*ca)
        R = np.array([[math.cos(0.5), 0, math.sin(0.5)], [0, 1, 0], [-math.sin(0.5), 0, math.cos(0.5)]]) @ \
            np.array([[1, 0, 0], [0, math.cos(0.8), -math.sin(0.8)], [0, math.sin(0.8), math.cos(0.8)]])
        d1 = R @ np.array([sa, 0, ca]); d2 = R @ np.array([-sa, 0, ca])
        z0 = np.concatenate([spinor(d1, s_len, 0.37), spinor(d2, s_len, -1.21)])
        s1, s2 = spins(z0); S = s1 + s2; n = S/np.linalg.norm(S); Sz = S[2]
        res = {}
        for lab, Mg, want in (('C_a', M_Na, 2*s_len), ('Z', M_Sz + 0.5*(M_Na + M_Nb), Sz + 2*s_len), ('G', M_Sn(n) + 0.5*(M_Na + M_Nb), np.linalg.norm(S) + 2*s_len)):
            mu, act, clos = cycle(z0, Mg, n); res[lab] = mu
            okA &= abs(act/(2*np.pi) - want) < 1e-9*max(1, want) and clos < 1e-9
            P('  |s| %5.2f |S| %6.3f  %-3s: action/2pi %.10f (predicted %.10f), closure %.1e, Maslov mu = %+.6f' % (s_len, np.linalg.norm(S), lab, act/(2*np.pi), want, clos, mu))
        muG.append(res['G']); muZ.append(res['Z'])
P('H730a (actions = 2 pi x predicted to 1e-9, cycles closed): %s' % okA)
mG = sorted(set(round(m) % 4 for m in muG))
P('H730b (report): mu(C_a) = 4 expected; mu(G) mod 4 over all points = %s; mu(Z) mod 4 = %s' % (mG, sorted(set(round(m) % 4 for m in muZ))))

# ---- H730c: spectrum for spins 10 + 10
P('\nH730c: s1 = s2 = 10 (|s| = 10.5): EBK |S| from mu(G) vs QM; E/A = (|S|^2 - 2 |s|^2)/2 vs [S(S+1) - 2 s(s+1)]/2')
s = 10; sl = s + 0.5; QM = [(S*(S + 1) - 2*s*(s + 1))/2 for S in range(0, 2*s + 1)]
for lab, off in (('|S| = S + 1/2 (mu(G) = 2 mod 4)', 0.5), ('|S| = S + 1   (mu(G) = 0 mod 4)', 1.0)):
    E = [((S + off)**2 - 2*sl*sl)/2 for S in range(0, 2*s + 1)]; err = [e - q for e, q in zip(E, QM)]
    P('  %s: error min %+.4f max %+.4f (at S = 0: %+.4f, S = 20: %+.4f)' % (lab, min(err), max(err), err[0], err[-1]))
m = mG[0] if len(mG) == 1 else None
if m is not None:
    allowed = [x for x in np.arange(0.25, 2*sl, 0.25) if abs((x + 2*sl - m/4) - round(x + 2*sl - m/4)) < 1e-9]
    P('  EBK from the computed mu(G) = %d (mod 4): allowed |S| in (0, %.1f): %s' % (m, 2*sl, ', '.join('%.2f' % x for x in allowed)))
    okC = (m == 2)
else:
    okC = False
P('H730c (EBK lengths from mu(G) give QM spectrum to O(hbar^2), i.e. |S| = S + 1/2): %s' % okC)
open(os.path.join(HERE, 'out', 'ana730.txt'), 'w').write('\n'.join(out) + '\n')

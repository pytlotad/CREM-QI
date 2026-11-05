# Audit 715 (offline): the zitterbewegung clock (2 m c^2/hbar) as the source of a FREE field with threshold 2 m c^2 -- memory response on Kepler
# orbits (hydrogen, infinite nucleus, atomic units) against the Compton clock / threshold m c^2 control (712 W1); analytic parts: hydrogen with a
# proton (electron-only vs joint clock against the measured 1S-2S) and positronium (single-particle vs joint clock).
# usage: python3 ana715.py
import os, math
import numpy as np
from multiprocessing import Pool

M = 128; K = 400
NU = np.round(np.arange(0.80, 4.6001, 0.01), 2)
VAR = {'C1': (1, 1.0), 'Z2': (2, 2.0)}      # (q: clock multiple of m c^2, M_f: field mass in m)
out = []


def kepler(nu, L):
    a = nu*nu; e = math.sqrt(max(0.0, 1 - (L/nu)**2)); T = 2*math.pi*nu**3
    Mean = (np.arange(M))*2*math.pi/M; u = Mean.copy()
    for _ in range(60): u = u - (u - e*np.sin(u) - Mean)/(1 - e*np.cos(u))
    x = a*(np.cos(u) - e); y = a*math.sqrt(1 - e*e)*np.sin(u)
    return x, y, u, e, T


def response(args):
    W, nu, L = args
    q, Mf = VAR[W]
    x, y, u, e, T = kepler(nu, L); dt = T/M
    thp = -q*nu*e*np.sin(u); D0 = -q/(2*nu*nu)
    R = 0j
    idx = np.arange(M)
    for s in range(M):                                   # shift s: tau = (s + m M) dt
        j = (idx - s) % M
        dx2 = (x - x[j])**2 + (y - y[j])**2; dth = thp - thp[j]
        m0 = 1 if s < M//2 else 0                         # tau > T/2 only (memory part)
        tau = (s + M*np.arange(m0, K))*dt
        if tau.size == 0: continue
        ph = np.exp(1j*(Mf*dx2[None, :]/(2*tau[:, None]) + dth[None, :])).mean(1)
        R += np.sum(tau**-1.5*np.exp(1j*D0*tau)*ph)*dt
    return W, nu, L, R*math.sqrt(T)                       # normalised to returns: tau^{-3/2} dtau = T^{-1/2} k^{-3/2} dk


def maxima(nus, A, lo):
    return [nus[i] for i in range(1, len(A) - 1) if A[i] > A[i - 1] and A[i] >= A[i + 1] and nus[i] > lo]


def main():
    jobs = [(W, nu, L) for W in VAR for L in (0.5, 1.0, 1.5) for nu in NU if nu > L + 0.005]
    with Pool(4) as p: res = p.map(response, jobs)
    tab = {}
    for W, nu, L, R in res: tab.setdefault((W, L), []).append((nu, R))
    okA = True; okB = True; strength2 = {}
    out.append('PART 1: memory response |R| (normalised to returns, sum k^-3/2 -> zeta(3/2) = 2.612 for a pure tone at threshold), tau > T/2, 400 periods')
    for (W, L), rows in sorted(tab.items()):
        nus = np.array([r[0] for r in rows]); A = np.abs(np.array([r[1] for r in rows]))
        mx = maxima(nus, A, L + 0.02)
        step = 1 if W == 'Z2' else 2
        expect = [j for j in range(step, 5, step) if j > L + 0.02]
        near = all(min(abs(m - j) for j in range(1, 5)) < 0.03 for m in mx) if W == 'Z2' else all(min(abs(m - j) for j in (2, 4)) < 0.03 for m in mx)
        have = all(any(abs(m - j) < 0.03 for m in mx) for j in expect)
        if W == 'Z2': okA &= near and have
        else: okA &= near and have
        line = ', '.join(f'{m:.2f} (|R| {A[list(nus).index(m)]:.3f})' for m in mx)
        out.append(f'  {W} L {L}: local maxima at nu = {line}; expected {expect}: positions {near}, all present {have}')
        for j0, j1 in zip(expect[:-1], expect[1:]):
            i0 = int(np.argmin(abs(nus - j0))); i1 = int(np.argmin(abs(nus - j1))); mn = A[i0:i1 + 1].min()
            c = A[i0]/mn; okB &= c > 3.41
            out.append(f'     contrast |R(nu={j0})| / min between {j0} and {j1}: {c:.3f}')
        if W == 'Z2': strength2[L] = A[int(np.argmin(abs(nus - 2)))]
        k = int(np.argmin(abs(nus - 2.5)))
        out.append(f'     |R| at nu = 1.5, 2.5, 3.5 (between): ' + ', '.join(f'{A[int(np.argmin(abs(nus - v)))]:.3f}' for v in (1.5, 2.5, 3.5) if v > L))
    okC = strength2[0.5] > strength2[1.0] > strength2[1.5]
    out.append(f'H715a (Z2 maxima at nu = j, all j = 1-4 present; C1 at nu = 2, 4 only): {okA}')
    out.append(f'H715b (contrast > 3.41 between neighbouring resonances): {okB}')
    out.append(f'H715c (Z2 |R(nu = 2)|: L 0.5 {strength2[0.5]:.3f} > L 1.0 {strength2[1.0]:.3f} > L 1.5 {strength2[1.5]:.3f}): {okC}')
    # PART 2: analytic, hydrogen with a proton and positronium
    HZ = 6.579683920502e15; HEV = 27.211386245988; MP = 1836.15267343; C = 299792458.0
    def three(ha):
        f = ha*HZ; return f'{f/1e6:,.3f} MHz = {ha*HEV:.6f} eV = {C/f*1e9:.4f} nm'
    mu = 1/(1 + 1/MP); bohr = 0.75*0.5*mu; meas = 2466061413.187e6/HZ
    el = bohr/(1 + 1/MP)**2
    out.append('PART 2 (analytic): resonance of the free threshold-2mc^2 field: (deficit of the clock) + j hbar Omega = 0, Omega = 2 Ry_mu/(hbar nu^3)')
    out.append(f'  hydrogen, joint clock (deficit 2E exactly, 712): nu = j -> 1S-2S {three(bohr)} (Bohr, reduced mass)')
    out.append(f'  hydrogen, electron-only clock (deficit 2 (m_p/M) E): nu = j (1 + m_e/m_p) -> 1S-2S {three(el)}')
    out.append(f'  measured 1S-2S {three(meas)} (Parthey et al. 2011): measured - Bohr {(meas - bohr)*HZ/1e6:+,.1f} MHz, electron-only - measured {(el - meas)*HZ/1e6:+,.1f} MHz'
               f' ({abs(el - meas)/abs(meas - bohr):.0f}x)')
    okD = abs(el - meas) > 100*abs(meas - bohr)
    out.append(f'H715d (electron-only clock off by > 100x (measured - Bohr)): {okD}')
    ps = 0.75*0.25; psm = 1233607216.4e6/HZ
    out.append(f'  Ps joint clock: levels -{0.25*HEV:.4f}/n^2 eV, 1S-2S {three(ps)}; single-particle clock (deficit 2 x E/2 = E): nu = 2j, levels'
               f' -{0.25*HEV/4:.4f}/j^2 eV, 1S-2S {three(ps/4)}; measured {three(psm)} (Fee et al. 1993)')
    out.append(f'  Ps per channel: the clock deficit and the threshold condition do not depend on spin -> p-Ps (1^1S0 - 2^1S0) and o-Ps (1^3S1 - 2^3S1)'
               ' the same at this order (measured line: o-Ps)')
    okE = abs(ps/4 - psm)/psm > 0.5 and abs(ps - psm)/psm < 1e-3
    out.append(f'H715e (Ps: single-particle clock factor 4 off, joint clock = Bohr within 1e-3 of measured): {okE}')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana715.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

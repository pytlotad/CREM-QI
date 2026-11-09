# Audit 764 (offline): second Schwarzian iteration (hbar^4) of the EBK wave with a smoothed eta_k -> <r^-3>, fine structure.  -> out/ana764.txt
# Loop = ana737.u_iter with one change: eta_k = 1/4 + r^2 D_k (k >= 1) smoothed by a degree-5 least-squares spline in x = ln r (see prereg764.txt).
import os, sys, io, contextlib, math, time, numpy as np
from scipy.interpolate import CubicSpline, make_lsq_spline
D0 = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(D0, '..', 'section-737'))
with contextlib.redirect_stdout(io.StringIO()): import ana737 as w
out = []; pr = lambda s: (print(s, flush=True), out.append(s))
ALPHA = 7.2973525643e-3; AE = 1.15965218059e-3; MP_ME = 1836.15267343; MU = 1 / (1 + 1 / MP_ME)
HZ = 6.579683920502e15; EV = 27.211386245988; c_cm = 2.99792458e10; MEAS = 9911.200 + 1057.845; XI_PS = 2741.682; ENG = 10968.704747
un = lambda f: f"{f:.3f} MHz = {f*1e6*EV/HZ*1e6:.4f} µeV = {c_cm/(f*1e6):.4f} cm"

def smoother(r, eta, N):
    x = np.log(r); t = np.r_[[x[0]] * 6, np.linspace(x[0], x[-1], N + 2)[1:-1], [x[-1]] * 6]
    s = make_lsq_spline(x, eta, t, k=5)
    return (lambda z: s(np.log(z))), np.max(np.abs(s(x) - eta))

def u_iter_s(n, l, r, k, N=None, smooth0=False, glob=False):
    E = -0.5 / n**2; nr = n - l - 1; L2 = (l + 0.5)**2
    QL = lambda x: 2 * E + 2 / x - L2 / x**2
    disc = math.sqrt(1 + 2 * E * L2); tp = ((1 - disc) / (-2 * E), (1 + disc) / (-2 * E))
    Qf, Q, Qp, Qpp = QL, QL(r), (-2 / r**2 + 2 * L2 / r**3), (4 / r**3 - 6 * L2 / r**4)
    sig, Gam, cc, ph, tp = w.mapping(Qf, r, nr, tp, glob); phs = [ph]; res = []
    for it in range(k):
        D, sp = w.schwarz(sig, Gam, cc, r, Q, Qp, Qpp, tp); eta = 0.25 + r * r * D
        if N is not None and (it >= 1 or smooth0):
            et, rmax = smoother(r, eta, N); res.append(rmax)
        else:
            et = CubicSpline(r, eta)
        Qf = (lambda QL_, e_: (lambda x: QL_(x) + e_(x) / x**2))(QL, et)
        Q = QL(r) + et(r) / r**2; Qp = Qpp = None
        sig, Gam, cc, ph, tp = w.mapping(Qf, r, nr, tp, glob); phs.append(ph)
    u = w.wave(sig, r, nr); uq = w.u_qm(n, l, r)
    if np.trapezoid(u * uq, r) < 0: u = -u
    return u, uq, phs, res

def fs(r3, g=2 * (1 + AE)): return 0.5 * ALPHA**2 * r3 * 1.5 * MU**3 * ((g - 1) + g / MP_ME) * HZ / 1e6
def xi(r3): return (3 / 16) * ALPHA**2 * r3 * ((2 * 2 * (1 + AE) - 1) / 3) * HZ / 1e6

if __name__ == '__main__':
    t0 = time.time(); ref762 = {(2, 1): -6.579e-05}
    pr('Audit 764: EBK wave, second Schwarzian iteration with smoothed eta_k (LSQ spline deg 5 in ln r, N interior knots); variant "in"')
    tab = {}
    for n, l in [(2, 1), (3, 1), (3, 2), (4, 1), (4, 2), (4, 3)]:
        r = w.grid(n); q = 1 / (n**3 * l * (l + 0.5) * (l + 1)); row = {}
        runs = [('k0', 0, None, False), ('k1', 1, None, False), ('k1+K1', 1, 600, True), ('k2-N600', 2, 600, False)]
        if (n, l) == (2, 1): runs += [('k2-N300', 2, 300, False), ('k2-N1200', 2, 1200, False)]
        for name, k, N, s0 in runs:
            u, uq, phs, res = u_iter_s(n, l, r, k, N, s0)
            v = np.trapezoid(u * u / r**3, r); ov = 1 - np.trapezoid(u * uq, r)
            row[name] = (v / q - 1, ov, phs, v)
            pr(f"  {n}{'spdf'[l]} {name:9s}: <r^-3>/QM-1 {v/q-1:+.4e}  1-<u|uQM> {ov:.2e}  faza k: {' '.join(f'{p-1:+.2e}' for p in phs)}"
               + (f"  max|eta - splajn| {', '.join(f'{x:.1e}' for x in res)}" if res else '') + f"  [{time.time()-t0:.0f} s]")
        tab[(n, l)] = (q, row)
    pr('')
    q, row = tab[(2, 1)]
    k2 = [row[x][0] for x in ('k2-N300', 'k2-N600', 'k2-N1200')]
    K1 = abs(row['k1+K1'][0] - row['k1'][0]); K2 = abs(row['k1'][0] - ref762[(2, 1)])
    a = max(k2) - min(k2) < 1e-6 and K1 < 1e-7
    pr(f"K1 (eta_0 wygladzone, k = 1): zmiana {K1:.1e} (< 1e-7);  K2 (k = 1 = 762 do zaokraglenia 4 cyfr): {row['k1'][0]:+.4e} wobec {ref762[(2,1)]:+.3e}")
    pr(f"H764a: 2p k = 2 przy N = 300/600/1200: {' '.join(f'{x:+.4e}' for x in k2)}; rozrzut {max(k2)-min(k2):.1e} -> {'TAK' if a else 'NIE'}")
    b = abs(row['k2-N600'][0]) < 3e-6 and all(abs(rw['k2-N600'][0]) < abs(rw['k1'][0]) for _, rw in tab.values())
    pr(f"H764b: 2p k = 2 {row['k2-N600'][0]:+.3e}; lepiej niz k = 1 w {sum(abs(rw['k2-N600'][0]) < abs(rw['k1'][0]) for _, rw in tab.values())}/6 -> {'TAK' if b else 'NIE'}")
    pr(f"pomiar 2P3/2 - 2P1/2: {un(MEAS)};  silnik domyslny (763): {ENG:.3f} MHz;  <r^-3> QM: {un(fs(q))}")
    for name in ('k0', 'k1', 'k2-N600'):
        f = fs(row[name][3])
        pr(f"  {name}: {un(f)}; - QM {f-fs(q):+.3f} MHz; - pomiar {f-MEAS:+.3f} MHz ({(f/MEAS-1):+.2e}); - silnik {f-ENG:+.3f} MHz")
        pr(f"    o-Ps 2^3P_J stala LS: {xi(row[name][3]):.3f} MHz wobec {XI_PS} (teoria QED) ({xi(row[name][3])/XI_PS-1:+.2e}); p-Ps 2^1P_1: bez LS")
    c = abs(fs(row['k2-N600'][3]) - fs(q)) < 0.05
    pr(f"H764c: {'TAK' if c else 'NIE'}")
    pr(f"runtime {time.time()-t0:.0f} s")
    open(os.path.join(D0, 'out', 'ana764.txt'), 'w').write('\n'.join(out) + '\n')

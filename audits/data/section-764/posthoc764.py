# Audit 764 POST HOC: k = 2 with eta_1 = eta_0 (closed form) + S(eta_1 - eta_0), S = LSQ spline (deg 5, ln r, N knots) fitted with
# windows |r - r_t| < W r_t and r < 1e-5 removed (finite-difference D_1 is garbage there).  Control K3: same recipe applied at k = 1 to
# eta_0 from finite differences (schwarz with Qp = None) must reproduce k = 1.  -> out/posthoc764.txt
import os, sys, math, time, numpy as np
from scipy.interpolate import CubicSpline, make_lsq_spline
D0 = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D0)
import ana764 as a; w = a.w
out = []; pr = lambda s: (print(s, flush=True), out.append(s))
def wfit(r, d, tp, N, W):
    keep = (r >= 1e-5) & np.all([np.abs(r - t) >= W * t for t in tp], axis=0)
    x = np.log(r[keep]); t = np.r_[[x[0]] * 6, np.linspace(x[0], x[-1], N + 2)[1:-1], [x[-1]] * 6]
    s = make_lsq_spline(x, d[keep], t, k=5); res = s(x) - d[keep]
    return (lambda z: s(np.log(np.maximum(z, r[keep][0])))), np.sqrt(np.mean(res**2)), np.abs(res).max()
def run(n, l, r, k, N=None, W=None, K3=False):
    E = -0.5 / n**2; nr = n - l - 1; L2 = (l + 0.5)**2
    QL = lambda x: 2 * E + 2 / x - L2 / x**2
    disc = math.sqrt(1 + 2 * E * L2); tp = ((1 - disc) / (-2 * E), (1 + disc) / (-2 * E))
    Q, Qp, Qpp = QL(r), (-2 / r**2 + 2 * L2 / r**3), (4 / r**3 - 6 * L2 / r**4)
    sig, Gam, cc, ph, tp = w.mapping(QL, r, nr, tp); phs = [ph]; info = []; e0 = None
    for it in range(k):
        if it == 0:
            D, _ = w.schwarz(sig, Gam, cc, r, Q, Qp, Qpp, tp); e0 = 0.25 + r * r * D; base = CubicSpline(r, e0)
            if K3:
                Dfd, _ = w.schwarz(sig, Gam, cc, r, Q, None, None, tp); s, rms, mx = wfit(r, 0.25 + r * r * Dfd - e0, tp, N, W); info.append((rms, mx))
                et = (lambda b, s_: (lambda z: b(z) + s_(z)))(base, s)
            else: et = base
        else:
            D, _ = w.schwarz(sig, Gam, cc, r, Q, None, None, tp); s, rms, mx = wfit(r, 0.25 + r * r * D - e0, tp, N, W); info.append((rms, mx))
            et = (lambda b, s_: (lambda z: b(z) + s_(z)))(base, s)
        Qf = (lambda e_: (lambda x: QL(x) + e_(x) / x**2))(et); Q = QL(r) + et(r) / r**2
        sig, Gam, cc, ph, tp = w.mapping(Qf, r, nr, tp); phs.append(ph)
    u = w.wave(sig, r, nr); uq = w.u_qm(n, l, r)
    if np.trapezoid(u * uq, r) < 0: u = -u
    q = 1 / (n**3 * l * (l + 0.5) * (l + 1)); v = np.trapezoid(u * u / r**3, r)
    return v / q - 1, 1 - np.trapezoid(u * uq, r), phs, info, v, q
if __name__ == '__main__':
    t0 = time.time(); pr('Audit 764 POST HOC: k = 2 with windowed fit of eta_1 - eta_0 (N knots in ln r, windows W r_t around turning points)')
    fmt = lambda name, R: pr(f"  {name}: <r^-3>/QM-1 {R[0]:+.4e}  1-<u|uQM> {R[1]:.2e}  faza: {' '.join(f'{p-1:+.2e}' for p in R[2])}  "
                             + ' '.join(f'fit rms {i[0]:.1e} max {i[1]:.1e}' for i in R[3]) + f"  [{time.time()-t0:.0f} s]")
    r = w.grid(2); k1 = a.u_iter_s(2, 1, r, 1); q = 1 / (8 * 1.5 * 2); k1v = np.trapezoid(k1[0]**2 / r**3, r) / q - 1
    pr(f"2p k = 1 (odniesienie): {k1v:+.4e}")
    K3 = {}
    for N in (40, 80, 160):
        for W in (0.03, 0.06, 0.10):
            try: R = run(2, 1, r, 1, N, W, K3=True); K3[(N, W)] = R[0]; fmt(f"2p K3 k=1 N{N:3d} W{W:.2f}", R)
            except Exception as e: pr(f"  2p K3 N{N} W{W}: blad {e}")
    pr(f"K3: max |k1(K3) - k1| = {max(abs(v - k1v) for v in K3.values()):.1e}")
    S = {}
    for N in (40, 80, 160):
        for W in (0.03, 0.06, 0.10):
            try: R = run(2, 1, r, 2, N, W); S[(N, W)] = R; fmt(f"2p k=2 N{N:3d} W{W:.2f}", R)
            except Exception as e: pr(f"  2p k=2 N{N} W{W}: blad {e}")
    v2 = [R[0] for R in S.values()]; pr(f"2p k = 2: zakres {min(v2):+.4e} ... {max(v2):+.4e}, rozrzut {max(v2)-min(v2):.1e}; mediana {np.median(v2):+.4e}")
    R = S[(80, 0.06)]; f = a.fs(R[4]); fq = a.fs(R[5])
    pr(f"2P3/2 - 2P1/2 (N 80, W 0,06): {a.un(f)}; - QM {f-fq:+.3f}; - pomiar {f-a.MEAS:+.3f} ({f/a.MEAS-1:+.2e}); - silnik {f-a.ENG:+.3f} MHz;"
       f" rozrzut po (N, W): {(max(v2)-min(v2))*fq:.3f} MHz")
    pr(f"  o-Ps stala LS: {a.xi(R[4]):.3f} MHz wobec {a.XI_PS} ({a.xi(R[4])/a.XI_PS-1:+.2e}); p-Ps bez LS")
    for n, l in [(3, 1), (3, 2), (4, 1), (4, 2), (4, 3)]:
        r = w.grid(n)
        for N, W in ((80, 0.06), (160, 0.03)):
            try: fmt(f"{n}{'spdf'[l]} k=2 N{N:3d} W{W:.2f}", run(n, l, r, 2, N, W))
            except Exception as e: pr(f"  {n}{'spdf'[l]} N{N} W{W}: blad {e}")
    pr(f"runtime {time.time()-t0:.0f} s")
    open(os.path.join(D0, 'out', 'posthoc764.txt'), 'w').write('\n'.join(out) + '\n')

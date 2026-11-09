# Audit 766 (offline): diagnosis of the K3 error for p states and of the systematic p residual at k = 2 -> out/ana766.txt
import os, sys, io, math, time, contextlib, numpy as np
from scipy.interpolate import CubicSpline, make_lsq_spline
D0 = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(D0, '..', 'section-764'))
with contextlib.redirect_stdout(io.StringIO()): import posthoc764 as ph
w = ph.w
out = []; pr = lambda s: (print(s, flush=True), out.append(s))
def wfit(r, d, tp, N, W, rmin):
    keep = (r >= rmin) & np.all([np.abs(r - t) >= W * t for t in tp], axis=0)
    x = np.log(r[keep]); t = np.r_[[x[0]] * 6, np.linspace(x[0], x[-1], N + 2)[1:-1], [x[-1]] * 6]
    s = make_lsq_spline(x, d[keep], t, k=5)
    return (lambda z: s(np.log(np.maximum(z, r[keep][0]))))
def run(n, l, k, N=80, W=0.06, K3=False, fine=1, rmin=1e-5, glob=False):
    r = w.grid(n, fine); E = -0.5 / n**2; nr = n - l - 1; L2 = (l + 0.5)**2
    QL = lambda x: 2 * E + 2 / x - L2 / x**2
    disc = math.sqrt(1 + 2 * E * L2); tp = ((1 - disc) / (-2 * E), (1 + disc) / (-2 * E))
    Q, Qp, Qpp = QL(r), (-2 / r**2 + 2 * L2 / r**3), (4 / r**3 - 6 * L2 / r**4)
    sig, Gam, cc, ph_, tp = w.mapping(QL, r, nr, tp, glob)
    for it in range(k):
        if it == 0:
            D, _ = w.schwarz(sig, Gam, cc, r, Q, Qp, Qpp, tp); e0 = 0.25 + r * r * D; base = CubicSpline(r, e0); et = base
            if K3:
                Dfd, _ = w.schwarz(sig, Gam, cc, r, Q, None, None, tp); s = wfit(r, 0.25 + r * r * Dfd - e0, tp, N, W, rmin)
                et = (lambda b, s_: (lambda z: b(z) + s_(z)))(base, s)
        else:
            D, _ = w.schwarz(sig, Gam, cc, r, Q, None, None, tp); s = wfit(r, 0.25 + r * r * D - e0, tp, N, W, rmin)
            et = (lambda b, s_: (lambda z: b(z) + s_(z)))(base, s)
        Qf = (lambda e_: (lambda x: QL(x) + e_(x) / x**2))(et); Q = QL(r) + et(r) / r**2
        sig, Gam, cc, ph_, tp = w.mapping(Qf, r, nr, tp, glob)
    u = w.wave(sig, r, nr); uq = w.u_qm(n, l, r)
    q = 1 / (n**3 * l * (l + 0.5) * (l + 1))
    return np.trapezoid(u * u / r**3, r) / q - 1
if __name__ == '__main__':
    t0 = time.time(); ref = {}
    for L in open(os.path.join(D0, '..', 'section-765', 'out', 'ana765.txt')):
        p = L.split()
        if len(p) > 8 and p[0].endswith(':') and p[1] == 'k1': ref[p[0][:-1]] = (float(p[2]), float(p[4]), float(p[6]))
    pr('Audit 766: K3 and the p residual vs grid step (fine), fit cutoff (rmin) and phase-scale variant (in / gl); (W, N) = (0,06; 80)')
    okA = okB = okC = True; okCd = True
    for n, l in [(2, 1), (5, 1), (8, 1), (3, 2), (8, 2)]:
        nm = f"{n}{'spdf'[l]}"; pp = l == 1; T = lambda: f"[{time.time()-t0:.0f} s]"
        k1_1 = run(n, l, 1); K3_1 = run(n, l, 1, K3=True); k2_1 = run(n, l, 2)
        pr(f"  {nm} fine 1: k1 {k1_1:+.4e}  |K3-k1| {abs(K3_1-k1_1):.1e}  k2 {k2_1:+.4e}  (765: k1 {ref[nm][0]:+.3e}, K3 {ref[nm][1]:.1e}, k2 {ref[nm][2]:+.2e}) {T()}")
        k1_2 = run(n, l, 1, fine=2); K3_2 = run(n, l, 1, K3=True, fine=2); k2_2 = run(n, l, 2, fine=2)
        pr(f"  {nm} fine 2: k1 {k1_2:+.4e}  |K3-k1| {abs(K3_2-k1_2):.1e}  k2 {k2_2:+.4e}  (k1 zmiana {k1_2-k1_1:+.1e}, k2 zmiana {k2_2-k2_1:+.1e}) {T()}")
        a6 = run(n, l, 2, rmin=1e-6); a4 = run(n, l, 2, rmin=1e-4)
        pr(f"  {nm} rmin 1e-6: k2 {a6:+.4e} ({a6-k2_1:+.1e});  rmin 1e-4: k2 {a4:+.4e} ({a4-k2_1:+.1e}) {T()}")
        g1 = run(n, l, 1, glob=True); g2 = run(n, l, 2, glob=True)
        pr(f"  {nm} gl: k1 {g1:+.4e} ({g1-k1_1:+.1e})  k2 {g2:+.4e} (gl - in {g2-k2_1:+.1e}) {T()}")
        okA &= abs(K3_2 - k1_2) < 1e-8 and (not (n >= 5 and pp) or abs(K3_2 - k1_2) * 3 <= abs(K3_1 - k1_1))
        if pp:
            okB &= max(abs(k2_2 - k2_1), abs(a6 - k2_1), abs(a4 - k2_1)) < 3e-7; okC &= abs(g2 - k2_1) > 1e-6
        else: okCd &= abs(g2 - k2_1) < 3e-7
    pr(f"H766a (K3 przy fine 2 < 1e-8, spadek >= 3x dla 5p, 8p): {'TAK' if okA else 'NIE'}")
    pr(f"H766b (resztka p niezalezna od siatki i odciecia do 3e-7): {'TAK' if okB else 'NIE'}")
    pr(f"H766c (|gl - in| > 1e-6 dla p; d < 3e-7): p {'TAK' if okC else 'NIE'}, d {'TAK' if okCd else 'NIE'}")
    pr(f"runtime {time.time()-t0:.0f} s")
    open(os.path.join(D0, 'out', 'ana766.txt'), 'w').write('\n'.join(out) + '\n')

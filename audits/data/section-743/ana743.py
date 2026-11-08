# Audit 743 (M): H743a (delta identity), H743b (steady response scans 2D/3D), H743c (one-period free return phase), convergence.
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from modes743 import torus, modes2d, modes3d, delta, wrap, response
from scipy.optimize import brentq, minimize_scalar
pr = lambda *a: print(*a, flush=True)
beta = 0.2

pr("=== H743a: delta wlasnego modu na torusie EBK (2D beta 0,2) ===")
worst = 0
for nr, m in [(0, 1), (1, 0), (0, 2), (1, 1), (2, 1)]:
    T = torus(nr + 0.5, float(m), beta)
    md = [x for x in modes2d(beta, 3, 3) if x['nr'] == nr and x['m'] == m][0]
    d = delta(md, T); dev = abs(abs(wrap(d)) - np.pi); worst = max(worst, dev)
    pr(f"  (n_r, m) = ({nr}, {m}): E_EBK - E_n = {T['E'] - md['E']:+.3e}, delta = {d:.12f} = {d/np.pi:.12f} pi, |wrap| - pi = {dev:.2e}")
pr(f"  max ||wrap(delta)| - pi| = {worst:.2e}  (kryterium 1e-10)")
pr("=== H743a: 3D beta 0,2 (torus Langera-EBK, mod m = l) ===")
w3 = 0
for nr, l in [(0, 1), (1, 1), (0, 2)]:
    T = torus(nr + 0.5, l + 0.5, beta)
    md = [x for x in modes3d(beta, 3, 3) if x['nr'] == nr and x['l'] == l and x['m'] == l][0]
    d = delta(md, T); pred = np.pi * (1 + (l + 0.5) / T['Lam'])
    w3 = max(w3, abs(wrap(d - pred)))
    pr(f"  (n_r, l) = ({nr}, {l}): E_EBK - E_n = {T['E'] - md['E']:+.3e}, delta = {d/np.pi:.12f} pi, przewidywane {pred/np.pi:.12f} pi, "
       f"niedostrojenie od 0 mod 2pi = {abs(wrap(d)):.4f} rad")
pr(f"  max |delta - pi(1 + (l+1/2)/Lambda)| = {w3:.2e}")
pr("=== H743a: Kepler 3D (beta 0) przy J = n ===")
for Jr, L, nr, l in [(0.5, 1.5, 0, 1), (1.5, 1.5, 1, 1), (0.5, 2.5, 0, 2)]:
    T = torus(Jr, L, 0.0); md = [x for x in modes3d(0.0, 3, 3) if x['nr'] == nr and x['l'] == l and x['m'] == l][0]
    pr(f"  J_r {Jr}, L {L} (J = {Jr+L}): delta = {wrap(delta(md, T)):+.3e} rad")
pr("  Kepler 2D (beta 0) przy torusie EBK J_r 0,5 L 1 (mod n_r 0, m 1):",
   f"{delta([x for x in modes2d(0.0, 2, 2) if x['nr']==0 and x['m']==1][0], torus(0.5, 1.0, 0.0))/np.pi:.12f} pi")

# ---------- convergence ----------
T0 = torus(0.5, 1.0, beta); G = 0.01 / T0['Tr']
pr("=== zbieznosc P (2D, torus (0,6; 1,0), Gamma T_r0 = 0,01) ===")
for (nrm, mm, M) in [(8, 4, 2048), (12, 6, 4096), (16, 8, 4096), (12, 6, 8192)]:
    P, Z = response(torus(0.6, 1.0, beta, M), modes2d(beta, nrm, mm), G)
    pr(f"  nrmax {nrm} mmax {mm} M {M}: P = {P:.8e}, |Z| = {abs(Z):.6e}")
M2 = modes2d(beta, 12, 6)

def scan(xs, f):
    return np.array([f(x) for x in xs])

def maxima(xs, ys, f):
    out = []
    for i in range(1, len(xs) - 1):
        if ys[i] > ys[i - 1] and ys[i] >= ys[i + 1]:
            res = minimize_scalar(lambda x: -f(x), bounds=(xs[i - 1], xs[i + 1]), method='bounded', options=dict(xatol=1e-7))
            out.append((res.x, -res.fun))
    return out

def roots(xs, mods, tor):
    """roots of delta_n = 0 mod 2 pi along the scan (sign change of wrap without the +-pi jump)."""
    rts = []
    for md in mods:
        if abs(md['m']) > 4: continue
        w = np.array([wrap(delta(md, tor(x))) for x in xs])
        for i in range(len(xs) - 1):
            if w[i] * w[i + 1] < 0 and abs(w[i] - w[i + 1]) < np.pi:
                x0 = brentq(lambda x: wrap(delta(md, tor(x))), xs[i], xs[i + 1], xtol=1e-12)
                rts.append((x0, md['nr'], md['m'], md.get('l')))
    return sorted(rts)

def report(name, xs, tor, mods, ebk):
    f = lambda x: response(tor(x), mods, G)[0]
    ys = scan(xs, f); mx = maxima(xs, ys, f); rts = roots(xs, mods, tor)
    Pe = f(ebk); i0 = np.argmin(abs(xs - ebk))
    loc = ys[i0] >= ys[max(i0 - 1, 0)] and ys[i0] >= ys[min(i0 + 1, len(xs) - 1)]
    pr(f"--- {name}: P(EBK = {ebk}) = {Pe:.6e}; EBK lokalnym maksimum: {'TAK' if loc else 'NIE'}; P max skanu {ys.max():.6e} przy {xs[ys.argmax()]:.4f}")
    allok = True
    for x, p in mx:
        near = min(rts, key=lambda r: abs(r[0] - x)) if rts else None
        dist = abs(near[0] - x) if near else np.inf
        allok &= dist < 0.01
        pr(f"    maksimum P przy {x:.6f} (P = {p:.6e}, P/P_EBK = {p/Pe:.3f}); od EBK {x - ebk:+.6f}; najblizszy pierwiastek delta_n = 0: "
           f"{near[0]:.6f} (n_r {near[1]}, m {near[2]}{'' if near[3] is None else ', l %d' % near[3]}), odleglosc {dist:.2e}")
    pr(f"    pierwiastki delta_n = 0 w skanie: " + ", ".join(f"{r[0]:.4f}({r[1]},{r[2]}{'' if r[3] is None else ','+str(r[3])})" for r in rts))
    return ys, mx, allok, loc

pr("=== H743b 2D (beta 0,2, Gamma T_r0 = 0,01, nrmax 12, mmax 6, M 4096) ===")
xs = np.round(np.arange(0.1, 1.0 + 1e-9, 0.0025), 6)
ysA, mxA, okA, locA = report("skan J_r przy L = 1", xs, lambda x: torus(x, 1.0, beta), M2, 0.5)
xs2 = np.round(np.arange(0.6, 1.4 + 1e-9, 0.0025), 6)
ysB, mxB, okB, locB = report("skan L przy J_r = 0,5", xs2, lambda x: torus(0.5, x, beta), M2, 1.0)
pr(f"H743b 2D: EBK nie maksimum (oba skany): {'TAK' if not locA and not locB else 'NIE'}; maksima w 0,01 od delta_n = 0: {'TAK' if okA and okB else 'NIE'}")
np.savetxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', 'scan2d_Jr.txt'), np.c_[xs, ysA], header='J_r P (L=1, beta 0.2)')
np.savetxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', 'scan2d_L.txt'), np.c_[xs2, ysB], header='L P (J_r=0.5, beta 0.2)')

pr("=== H743b 3D (beta 0,2, nrmax 10, lmax 6) ===")
M3 = modes3d(beta, 10, 6); T3 = torus(0.5, 1.5, beta); G = 0.01 / T3['Tr']
ys3A, mx3A, _, loc3A = report("3D skan J_r przy L = 1,5", xs, lambda x: torus(x, 1.5, beta), M3, 0.5)
xs3 = np.round(np.arange(1.1, 1.9 + 1e-9, 0.0025), 6)
ys3B, mx3B, _, loc3B = report("3D skan L przy J_r = 0,5", xs3, lambda x: torus(0.5, x, beta), M3, 1.5)
best = [m for m in mx3A + mx3B]
pr("H743b 3D: maksima w 1e-3 od EBK: " + ("TAK" if any(abs(x - 0.5) < 1e-3 for x, _ in mx3A) or any(abs(x - 1.5) < 1e-3 for x, _ in mx3B) else "NIE"))
np.savetxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', 'scan3d_Jr.txt'), np.c_[xs, ys3A], header='J_r P (3D, L=1.5, beta 0.2)')
np.savetxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', 'scan3d_L.txt'), np.c_[xs3, ys3B], header='L P (3D, J_r=0.5, beta 0.2)')

pr("=== H743c (M): zrodlo przez jeden okres, Gamma = 0; arg(Z_{k+1}/Z_k) ===")
for name, Jr, L, b in [("EBK (0,5; 1)", 0.5, 1.0, beta), ("(0,6; 1,1)", 0.6, 1.1, beta), ("(0,4; 0,9)", 0.4, 0.9, beta), ("Kepler (0,5; 1)", 0.5, 1.0, 0.0)]:
    Zk = response(torus(Jr, L, b), modes2d(b, 12, 6), 0.0, onep=True)
    a1 = np.angle(Zk[1] / Zk[0]); a2 = np.angle(Zk[2] / Zk[1])
    pr(f"  {name}: |Z_1..3| = {abs(Zk[0]):.4e} {abs(Zk[1]):.4e} {abs(Zk[2]):.4e}; arg(Z2/Z1) = {a1:+.4f} rad, arg(Z3/Z2) = {a2:+.4f} rad; "
       f"|arg| - pi = {abs(a1) - np.pi:+.4f}, {abs(a2) - np.pi:+.4f}")

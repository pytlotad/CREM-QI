# Audit 748 (1): field-equation spectrum by finite differences -> classical radial action J_r at E_n.
import numpy as np, sys
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import quad
from scipy.optimize import brentq
pr = lambda *a: print(*a, flush=True)
POT = {'P1': lambda r: -1 / r + 0.2 / r**2, 'P2': lambda r: -1 / r + 0.1 * r, 'P3': lambda r: 0.5 * r**2 + 0.1 * r**4}
RMAX = {'P1': 500.0, 'P2': 60.0, 'P3': 12.0}
def spectrum(V, cen, n, rmax, N):
    r = np.linspace(rmax / N, rmax, N); h = r[1] - r[0]
    d = 1 / h**2 + cen / (2 * r**2) + V(r); e = -0.5 / h**2 * np.ones(N - 1)
    return eigh_tridiagonal(d, e, select='i', select_range=(0, n - 1), eigvals_only=True)
def Jr(V, E, L):
    W = lambda r: 2 * (E - V(r)) - L * L / r**2
    rs = np.geomspace(1e-4, 200, 200000); w = W(rs); idx = np.where(np.diff(np.sign(w)) != 0)[0]
    r1 = brentq(W, rs[idx[0]], rs[idx[0] + 1], xtol=1e-14); r2 = brentq(W, rs[idx[-1]], rs[idx[-1] + 1], xtol=1e-14)
    # substitution r = r1 + (r2 - r1)(1 - cos t)/2 removes the square-root endpoints
    f = lambda t: np.sqrt(max(W(r1 + (r2 - r1) * (1 - np.cos(t)) / 2), 0)) * (r2 - r1) * np.sin(t) / 2
    return quad(f, 0, np.pi, epsabs=1e-13, epsrel=1e-13, limit=400)[0] / np.pi
def exactP1(nr, Lq, dim):
    nu = np.sqrt(Lq + 2 * 0.2); N = nr + nu + (0.5 if dim == 2 else 0.5)  # Lq = m^2 (2D) or (l+1/2)^2 (3D): N = n_r + 1/2 + sqrt(Lq + 2 beta)
    return -0.5 / N**2
NR = 7
res = {}
for P in ['P1', 'P2', 'P3']:
    V = POT[P]
    for dim, q in [(2, 0), (2, 1), (2, 2), (2, 3), (3, 0), (3, 1), (3, 2), (3, 3)]:
        if P != 'P1' and q != 1: continue
        cen = (q * q - 0.25) if dim == 2 else q * (q + 1)
        NP = 125000 if P == 'P1' else 30000
        E1 = spectrum(V, cen, NR, RMAX[P], NP); E2 = spectrum(V, cen, NR, RMAX[P], 2 * NP)
        Ls = [abs(q)] if dim == 2 else [q + 0.5, np.sqrt(q * (q + 1))]
        for L in Ls:
            tag = f"{P} {dim}D {'m' if dim == 2 else 'l'}={q} L={L:.4f}" + ('' if dim == 2 else (' (Langer)' if L == q + 0.5 else ' (sqrt(l(l+1)))'))
            rows = []
            for nr in range(NR):
                E = E2[nr] + (E2[nr] - E1[nr]) / 3  # Richardson (h^2)
                if P == 'P1' and dim == 3 and L == 0 and q == 0: pass
                J = Jr(V, E, L) if L > 0 or P != 'P1' or True else np.nan
                ex = exactP1(nr, q * q if dim == 2 else (q + 0.5)**2, dim) if P == 'P1' else np.nan
                rows.append((nr, E, J, J - nr - 0.5, abs(E2[nr] - E1[nr]), (E / ex - 1) if P == 'P1' else np.nan))
            res[tag] = rows
            pr(f"--- {tag}")
            for nr, E, J, dv, dE, rel in rows:
                pr(f"   n_r {nr}: E = {E:+.10f}  J_r = {J:.8f}  J_r - n_r - 1/2 = {dv:+.2e}  (siatka 30k/60k dE {dE:.1e}" +
                   (f", E/scisle - 1 {rel:+.1e})" if P == 'P1' else ")"))
pr("=== OCENA ===")
a = max(abs(r[3]) for k, v in res.items() if k.startswith('P1') and ('2D' in k or 'Langer' in k) for r in v[:5])
pr(f"H748a: max |J_r - n_r - 1/2| (P1, 2D m 0-3, 3D Langer l 0-3, n_r 0-4) = {a:.2e} -> {'TAK' if a < 1e-4 else 'NIE'}")
b = [min(abs(r[3]) for r in v) for k, v in res.items() if k.startswith('P1 3D l=1 ') and 'sqrt' in k or k.startswith('P1 3D l=2 ') and 'sqrt' in k]
pr(f"H748b: min |J_r - n_r - 1/2| bez Langera, l = 1, 2: {', '.join(f'{x:.3e}' for x in b)} -> {'TAK' if min(b) > 1e-2 else 'NIE'}")
ok = True
for k, v in res.items():
    if k.startswith('P2') or k.startswith('P3'):
        if 'sqrt' in k: continue
        dv = [abs(r[3]) for r in v]; mono = all(dv[i + 1] <= dv[i] + 1e-9 for i in range(len(dv) - 1)); small = max(dv[3:]) < 0.02
        far = min(abs(r[2] - r[0]) for r in v) > 0.45; ok &= mono and small and far
        pr(f"  {k}: |dev| " + ' '.join(f"{x:.2e}" for x in dv) + f"; maleje {'TAK' if mono else 'NIE'}, n_r>=3 < 0,02 {'TAK' if small else 'NIE'}, od n_r > 0,45 {'TAK' if far else 'NIE'}")
pr(f"H748c: {'TAK' if ok else 'NIE'}")

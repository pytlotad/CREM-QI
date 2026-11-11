# Audit 779: closed form of the per-period drive phase Phi_n(J, L) of mode n = (n_r, m) on a frozen torus (derivation in prereg779.txt).
import numpy as np
def wrap(x): return (np.asarray(x) + np.pi) % (2 * np.pi) - np.pi
def Nmode(nr, m, beta): return nr + 0.5 + np.sqrt(m * m + 2 * beta)
def Phi(J, L, nr, m, beta, caust):
    Lam = np.sqrt(L * L + 2 * beta); N = J + Lam; Nn = Nmode(nr, m, beta)
    return wrap(2 * np.pi * J + 2 * caust + 2 * np.pi * L * (L - m) / Lam + np.pi * N - np.pi * N ** 3 / Nn ** 2)
def fejer(x, K):
    s = np.sin(np.asarray(x) / 2); out = np.full(np.shape(s), float(K * K)); nz = np.abs(s) > 1e-12
    out[nz] = np.sin(K * np.asarray(x)[nz] / 2) ** 2 / s[nz] ** 2; return out
def gammaT(GT, beta, N): N0 = 0.5 + np.sqrt(1 + 2 * beta); return GT * N ** 3 / N0 ** 3     # Gamma T = GT * T / Tr0

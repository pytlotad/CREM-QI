# Audit 755: Bohmian trajectories with a detector (pointer x, H_D = kappa p_x sigma_gg) for a two-level atom in N quantized modes (RWA, one excitation).
# Psi(q, x) = A_e(x,t)|e>Phi0(q) + |g>Phi0(q) sum_k B_k(x,t) sqrt(2 w_k) q_k.  Schrodinger on an x grid (FFT shift of the g branch, local coupling, Strang).
# Ensemble F (beables q and x):  dx/dt = kappa |S|^2/den,  dq_k/dt = [Im(S* u_k) - b_k Im(A_e* S)]/den,  den = |A_e|^2 + |S|^2  (at the actual x).
# Ensemble P (beable x only):     dx/dt = kappa rho_g(x)/(rho_e(x) + rho_g(x)).
# usage: python3 bohm755.py h M_F M_P out.npz
import sys, os, numpy as np
h, MF, MP, out = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
G = 0.01; W = 40 * G; w0 = 1.0; tau = 1 / G; T = float(os.environ.get("T755_TAU", "5")) * tau;  # T755_TAU only for quick tests
N = 200; kap = 0.15; Dl = 2.0
om = w0 - W / 2 + W * (np.arange(N) + 0.5) / N; dom = W / N; g = np.sqrt(G * dom / (2 * np.pi)); det = om - w0
sq = np.sqrt(2 * om); b = 2 * g / sq
Nx = 2048; x0g = -10.0; L = 96.0; dx = L / Nx; X = x0g + dx * np.arange(Nx); kx = 2 * np.pi * np.fft.fftfreq(Nx, dx)
shift = np.exp(-1j * kx * kap * h / 2)                     # half-step translation of the g branch by kappa h/2
Ae = (2 * np.pi) ** -0.25 * np.exp(-X ** 2 / 4) + 0j; B = np.zeros((N, Nx), complex)   # interaction picture (phases e^{-i w t} removed)
def coup(t, a, Bm):
    ph = np.exp(-1j * det * t); return -1j * g * (ph @ Bm), -1j * g * np.conj(ph)[:, None] * a[None, :]
def step(t, a, Bm):
    Bm = np.fft.ifft(np.fft.fft(Bm, axis=1) * shift[None, :], axis=1)
    k1 = coup(t, a, Bm); k2 = coup(t + h / 2, a + h / 2 * k1[0], Bm + h / 2 * k1[1]); k3 = coup(t + h / 2, a + h / 2 * k2[0], Bm + h / 2 * k2[1])
    k4 = coup(t + h, a + h * k3[0], Bm + h * k3[1])
    a = a + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); Bm = Bm + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
    return a, np.fft.ifft(np.fft.fft(Bm, axis=1) * shift[None, :], axis=1)
def interp(arr, x):   # linear interpolation along the last axis at positions x
    j = np.clip((x - x0g) / dx, 0, Nx - 1.000001); i = j.astype(int); f = j - i
    return arr[..., i] * (1 - f) + arr[..., i + 1] * f
xF = np.empty(MF); qF = np.empty((MF, N)); xP = np.empty(MP)
for m in range(max(MF, MP)):
    r = np.random.default_rng(m); x = r.normal(0, 1)
    if m < MP: xP[m] = x
    if m < MF: xF[m] = x; qF[m] = r.normal(0, 1, N) / np.sqrt(2 * om)
xF0, xP0 = xF.copy(), xP.copy()
def stateF(t, a, Bm, x, q):
    A = interp(a, x) * np.exp(-1j * w0 * t); u = (interp(Bm, x) * (np.exp(-1j * om * t) * sq)[:, None]).T     # (M,), (M, N)
    S = np.einsum('mk,mk->m', u, q); den = np.abs(A) ** 2 + np.abs(S) ** 2 + 1e-300
    vq = (np.imag(np.conj(S)[:, None] * u) - b[None, :] * np.imag(np.conj(A) * S)[:, None]) / den[:, None]
    w = np.abs(S) ** 2 / den; return kap * w, vq, w
def velP(re, rg, x):
    e, gg = interp(re, x), interp(rg, x); return kap * gg / (e + gg + 1e-300)
grid = np.arange(0, 5 + 1e-9, 0.05); gi = 1
WF = np.zeros((len(grid), MF), np.float32); XF = np.zeros((len(grid), MF), np.float32); XP = np.zeros((len(grid), MP), np.float32)
Pe = np.zeros(len(grid)); rhoS = {}
XF[0] = xF; XP[0] = xP; Pe[0] = 1.0; WF[0] = 0
trF = np.full(MF, -1.0); trP = np.full(MP, -1.0); nst = int(round(T / h)); t = 0.0
for s in range(nst):
    a1, B1 = step(t, Ae, B)
    am, Bmid = 0.5 * (Ae + a1), 0.5 * (B + B1)                            # midpoint (interaction picture, slow)
    re0, rg0 = np.abs(Ae) ** 2, (np.abs(B) ** 2).sum(0); re1, rg1 = np.abs(a1) ** 2, (np.abs(B1) ** 2).sum(0); rem, rgm = 0.5 * (re0 + re1), 0.5 * (rg0 + rg1)
    # P: RK4 in x
    k1 = velP(re0, rg0, xP); k2 = velP(rem, rgm, xP + h / 2 * k1); k3 = velP(rem, rgm, xP + h / 2 * k2); k4 = velP(re1, rg1, xP + h * k3)
    xPn = xP + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    # F: RK4 in (x, q)
    c1 = stateF(t, Ae, B, xF, qF); c2 = stateF(t + h / 2, am, Bmid, xF + h / 2 * c1[0], qF + h / 2 * c1[1])
    c3 = stateF(t + h / 2, am, Bmid, xF + h / 2 * c2[0], qF + h / 2 * c2[1]); c4 = stateF(t + h, a1, B1, xF + h * c3[0], qF + h * c3[1])
    xFn = xF + h / 6 * (c1[0] + 2 * c2[0] + 2 * c3[0] + c4[0]); qF = qF + h / 6 * (c1[1] + 2 * c2[1] + 2 * c3[1] + c4[1])
    for xo, xn, x0, tr in [(xP, xPn, xP0, trP), (xF, xFn, xF0, trF)]:
        new = (xn - x0 >= Dl) & (xo - x0 < Dl) & (tr < 0); tr[new] = (t + h * (Dl - (xo[new] - x0[new])) / (xn[new] - xo[new])) / tau
    xP, xF, Ae, B, t = xPn, xFn, a1, B1, t + h
    if gi < len(grid) and t >= grid[gi] * tau - 1e-9:
        XP[gi] = xP; XF[gi] = xF; WF[gi] = stateF(t, Ae, B, xF, qF)[2]; Pe[gi] = (np.abs(Ae) ** 2).sum() * dx
        if abs(grid[gi] - round(grid[gi])) < 1e-9: rhoS[int(round(grid[gi]))] = (np.abs(Ae) ** 2 + (np.abs(B) ** 2).sum(0)).astype(np.float64)
        gi += 1
np.savez_compressed(out, t=grid, X=X, Pe=Pe, XP=XP, XF=XF, WF=WF, trP=trP, trF=trF, xP0=xP0, xF0=xF0,
                    rho1=rhoS.get(1), rho2=rhoS.get(2), norm_end=(np.abs(Ae) ** 2 + (np.abs(B) ** 2).sum(0)).sum() * dx, h=h)
print(f"{out}: h {h} MF {MF} MP {MP}: P_e(5 tau) {Pe[-1]:.6e}, norma {(np.abs(Ae)**2 + (np.abs(B)**2).sum(0)).sum()*dx:.8f}", flush=True)

# Audit 754: Bohmian field trajectories for a two-level atom in N quantized field modes (RWA, single excitation; omega_0 = 1).
# Psi(q,t) = c_e |e> Phi0(q) + |g> S(q,t) Phi0(q),  S = sum_k u_k q_k,  u_k = c_k sqrt(2 omega_k).
# Guidance (field configuration q is the beable): dq_k/dt = [Im(S* u_k) - 2 g/sqrt(2 omega_k) Im(c_e* S)] / (|c_e|^2 + |S|^2 + |c_g0|^2);
# the second term is the interaction current of the RWA coupling, which contains d/dq_k in the q representation (prereg754, uzupelnienie 1).
# Start 'g' = |g, vacuum> (zero-excitation sector, amplitude c_g0 = 1, decoupled in RWA): adds |c_g0|^2 to den and to |g> weight, no velocity.
# usage: python3 bohm754.py N M start(e|g) dt_max out.npz [f_ad=0.5]
import sys, numpy as np
N, M, start, dt, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], float(sys.argv[4]), sys.argv[5]
G = 0.01; W = 40 * G; w0 = 1.0; tau = 1 / G; T = 5 * tau
om = w0 - W / 2 + W * (np.arange(N) + 0.5) / N; dom = W / N; g = np.sqrt(G * dom / (2 * np.pi))
f_ad = float(sys.argv[6]) if len(sys.argv) > 6 else 0.5   # adaptive step: dt = min(dt_max, f_ad |c_e|^2) (prereg754, uzupelnienie 2)
det = om - w0; sq = np.sqrt(2 * om); pg0 = 1.0 if start == 'g' else 0.0; b = 2 * g / sq
def rhs(t, ae, ak, q):
    ph = np.exp(-1j * det * t); dae = -1j * g * np.sum(ak * ph); dak = -1j * g * ae * np.conj(ph)
    u = ak * np.exp(-1j * om * t) * sq; ce = ae * np.exp(-1j * w0 * t); S = q @ u; den = abs(ce) ** 2 + np.abs(S) ** 2 + pg0
    v = (np.imag(np.conj(S)[:, None] * u[None, :]) - np.imag(np.conj(ce) * S)[:, None] * b[None, :]) / den[:, None]
    return dae, dak, v, (np.abs(S) ** 2 + pg0) / den
ae = 1.0 + 0j if start == 'e' else 0j; ak = np.zeros(N, complex)
q = np.empty((M, N))
for m in range(M): q[m] = np.random.default_rng(m).normal(0, 1, N) / np.sqrt(2 * om)
grid = np.arange(0, 5 + 1e-9, 0.05); ts = grid; Wt = np.zeros((len(grid), M)); Pe = np.zeros(len(grid))
te = np.full(M, -1.0); cross = np.zeros(M, int); vmax = 0.0; t = 0.0; nst = 0
k1 = rhs(t, ae, ak, q); wprev = k1[3]; Wt[0] = wprev; Pe[0] = 1.0 if start == 'e' else 0.0; gi = 1
while t < T - 1e-12:
    h = min(dt, f_ad * max(abs(ae) ** 2, 1e-300) if start == 'e' else dt, T - t)
    k1 = rhs(t, ae, ak, q); vmax = max(vmax, np.abs(k1[2]).max())
    k2 = rhs(t + h / 2, ae + h / 2 * k1[0], ak + h / 2 * k1[1], q + h / 2 * k1[2])
    k3 = rhs(t + h / 2, ae + h / 2 * k2[0], ak + h / 2 * k2[1], q + h / 2 * k2[2])
    k4 = rhs(t + h, ae + h * k3[0], ak + h * k3[1], q + h * k3[2])
    ae = ae + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); ak = ak + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
    q = q + h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]); t += h; nst += 1
    wn = rhs(t, ae, ak, q)[3]; pa, pb = wprev >= 0.5, wn >= 0.5; cross += pa != pb
    new = pb & ~pa & (te < 0); te[new] = (t - h + h * (0.5 - wprev[new]) / (wn[new] - wprev[new])) / tau
    while gi < len(grid) and grid[gi] * tau <= t + 1e-9:
        x = (grid[gi] * tau - (t - h)) / h; Wt[gi] = wprev + x * (wn - wprev); Pe[gi] = abs(ae) ** 2; gi += 1   # Pe at step end (h << tau)
    wprev = wn
np.savez_compressed(out, t=ts, w=Wt.astype(np.float32), Pe=Pe, te=te, cross=cross, vmax=vmax, N=N, M=M, dt=dt, f_ad=f_ad, steps=nst, start=start)
print(f"{out}: N {N} M {M} start {start} dt_max {dt} f {f_ad} steps {nst}: |c_e|^2(5 tau) {abs(ae)**2:.6e}, max|v| {vmax:.3e}", flush=True)

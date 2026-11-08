# Audit 743: steady / one-period response of the exact (bound-mode) matter-wave field driven by a point source on a prescribed
# classical torus of V = -1/r + beta/r^2, source phase phi_s = S0 - E t (gauge-coupled, 691/603) or S0/2 (rigid 446-like, PHASE K).
# hbar = mu = k = 1.  2D modes R(r) e^{i m phi}/sqrt(2 pi); 3D modes R_{n_r l}(r) Y_lm(pi/2, phi) (planar orbit in the equator).
import numpy as np
from scipy.special import eval_genlaguerre, gammaln, sph_harm_y

def torus(Jr, L, beta, M=4096):
    Lam = np.sqrt(L * L + 2 * beta); N = Jr + Lam; E = -0.5 / N**2; a = N * N; e = np.sqrt(max(0.0, 1 - Lam * Lam / (N * N)))
    Tr = 2 * np.pi * N**3; t = np.arange(M) * Tr / M; Mt = t / N**3
    u = Mt.copy()
    for _ in range(60): u -= (u - e * np.sin(u) - Mt) / (1 - e * np.cos(u))
    r = a * (1 - e * np.cos(u))
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(u / 2), np.sqrt(1 - e) * np.cos(u / 2))
    nu = np.unwrap(nu); nu[nu < 0] += 2 * np.pi * 0  # starts at 0 (periapsis)
    phi = (L / Lam) * nu; S0 = N * (u + e * np.sin(u)) - 2 * beta * nu / Lam
    dphi = 2 * np.pi * L / Lam; dS0 = 2 * np.pi * N - 4 * np.pi * beta / Lam
    return dict(Jr=Jr, L=L, beta=beta, Lam=Lam, N=N, E=E, a=a, e=e, Tr=Tr, t=t, r=r, phi=phi, S0=S0, dphi=dphi, dS0=dS0, M=M)

def modes2d(beta, nrmax, mmax):
    out = []
    for m in range(-mmax, mmax + 1):
        nu = np.sqrt(m * m + 2 * beta)
        for nr in range(nrmax + 1):
            N = nr + nu + 0.5; al = 2 * nu
            lognorm = 0.5 * (-2 * np.log(N / 2) - (gammaln(nr + al + 1) - gammaln(nr + 1)) - np.log(2 * N))
            out.append(dict(nr=nr, m=m, l=None, N=N, E=-0.5 / N**2, nu=nu, al=al, lognorm=lognorm, ang=1 / np.sqrt(2 * np.pi), dim=2))
    return out

def modes3d(beta, nrmax, lmax):
    out = []
    for l in range(lmax + 1):
        lam = -0.5 + np.sqrt((l + 0.5)**2 + 2 * beta)
        for m in range(-l, l + 1):
            if (l + m) % 2: continue
            Y = abs(complex(sph_harm_y(l, m, np.pi / 2, 0.0)))
            for nr in range(nrmax + 1):
                N = nr + lam + 1; al = 2 * lam + 1
                lognorm = 0.5 * (-3 * np.log(N / 2) - (gammaln(nr + al + 1) - gammaln(nr + 1)) - np.log(2 * N))
                out.append(dict(nr=nr, m=m, l=l, N=N, E=-0.5 / N**2, nu=lam, al=al, lognorm=lognorm, ang=Y, dim=3))
    return out

def radial(md, r):
    rho = 2 * r / md['N']
    return np.exp(md['lognorm'] + md['nu'] * np.log(rho) - rho / 2) * eval_genlaguerre(md['nr'], md['al'], rho)

def delta(md, T):
    return T['dS0'] - md['m'] * T['dphi'] + (md['E'] - T['E']) * T['Tr']

def wrap(x): return (x + np.pi) % (2 * np.pi) - np.pi

def response(T, mods, Gamma, phase='L', A=1.0, onep=False):
    """Steady state (onep=False): P = <sum|c_n|^2>, Z = <exp(-i phi_s) psi(x_p)>.  onep=True: source on for one period, Gamma = 0,
    returns Z_k (k = 1..3) for the free field."""
    t, Tr = T['t'], T['Tr']
    phis = (T['S0'] - T['E'] * t) if phase == 'L' else 0.5 * T['S0']
    dphis = (T['dS0'] - T['E'] * Tr) if phase == 'L' else 0.5 * T['dS0']
    dt = Tr / T['M']; P = 0.0; Z = 0j; Zk = np.zeros(3, complex)
    for md in mods:
        R = radial(md, T['r']) * md['ang']
        psi = R * np.exp(1j * md['m'] * T['phi'])           # psi_n(x_p(t)) without exp(-i E_n t)
        f = np.conj(psi) * np.exp(1j * (phis + md['E'] * t))
        d = dphis - md['m'] * T['dphi'] + md['E'] * Tr        # f(t + Tr) = f(t) exp(i d)
        if onep:
            G = np.sum(f) * dt; c = -1j * A * G
            h = np.sum(np.exp(-1j * phis) * psi * np.exp(-1j * md['E'] * t)) * dt / Tr
            for k in range(1, 4): Zk[k - 1] += c * h * np.exp(-1j * k * d)  # Z_k: factor exp(-i k d) per period
            continue
        w = np.exp(Gamma * t) * f
        G = np.concatenate([[0], np.cumsum(0.5 * (w[1:] + w[:-1]) * dt)])
        GT = G[-1] + 0.5 * (w[-1] + np.exp(Gamma * Tr) * f[0] * np.exp(1j * d)) * dt
        q = np.exp(-Gamma * Tr - 1j * d)
        I = np.exp(-Gamma * t) * (G + q * (GT - G))
        c = -1j * A * I / (1 - q)
        P += np.mean(np.abs(c)**2)
        Z += np.mean(np.exp(-1j * phis) * c * psi * np.exp(-1j * md['E'] * t))
    return (Zk if onep else (P, Z))

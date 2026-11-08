# Audit 751 H751d: populations driven by the gauge-phase source (743 steady state) on the Langer 2p torus (Kepler 3D, J_r 0.5, L 1.5).
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-743'))
from modes743 import torus, modes3d, radial
T = torus(0.5, 1.5, 0.0); G = 0.1 / T['Tr']; t = T['t']; Tr = T['Tr']; dt = Tr / T['M']
phis = T['S0'] - T['E'] * t; dphis = T['dS0'] - T['E'] * Tr
pop = {}
for md in modes3d(0.0, 10, 6):
    psi = radial(md, T['r']) * md['ang'] * np.exp(1j * md['m'] * T['phi'])
    f = np.conj(psi) * np.exp(1j * (phis + md['E'] * t)); d = dphis - md['m'] * T['dphi'] + md['E'] * Tr
    w = np.exp(G * t) * f; Gc = np.concatenate([[0], np.cumsum(0.5 * (w[1:] + w[:-1]) * dt)])
    GT = Gc[-1] + 0.5 * (w[-1] + np.exp(G * Tr) * f[0] * np.exp(1j * d)) * dt; q = np.exp(-G * Tr - 1j * d)
    I = np.exp(-G * t) * (Gc + q * (GT - Gc)); c = -1j * I / (1 - q)
    pop[(md['nr'], md['l'], md['m'])] = np.mean(np.abs(c)**2)
tot = sum(pop.values()); p1s = pop[(0, 0, 0)]; p2p = pop[(0, 1, 1)] + pop[(0, 1, -1)]; p2s = pop[(1, 0, 0)]
print(f"torus Langera 2p (J_r 0,5, L 1,5, e {T['e']:.4f}), Gamma T_r = 0,1, mody 3D n_r <= 10, l <= 6")
print(f"udzialy w sumie |c|^2: 2p(m=+1) {pop[(0,1,1)]/tot:.4f}, 2p(m=-1) {pop[(0,1,-1)]/tot:.4f}, 1s {p1s/tot:.4f}, 2s {p2s/tot:.4f}")
for k, v in sorted(pop.items(), key=lambda kv: -kv[1])[:8]: print(f"  (n_r {k[0]}, l {k[1]}, m {k[2]:+d}) {v/tot:.4f}")
s0 = p1s / (p1s + pop[(0, 1, 1)]); print(f"s0 = P(1s)/(P(1s) + P(2p, m=+1)) = {s0:.6e}")
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', 's0.txt'), 'w').write(f"{s0:.6e}\n")

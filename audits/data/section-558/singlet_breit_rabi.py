"""Ps 1S w polu: drugi punkt stały energii z audytu 539 to stan S' Breita–Rabiego.

Bez silnika, tylko numpy. Energia klasycznej pary spinów (n1 — e-, n2 — e+,
wektory jednostkowe, B wzdłuż z) z regułą 539 (Zeeman pary ważony sqrt(P_T)):
    E = K n1.n2 + Z (|n1 + n2| / 2) (n1 - n2).z,   K = dE_HFS/4,  x = Z/K.
QM w sektorze m = 0 ({S, T0}): H = [[-3K, 2Z], [2Z, K]].
Waga singletu konfiguracji klasycznej: w = (1 - n1.n2)/2 (w = 1 antyrównoległe).
Sprawdzamy: (1) oba punkty stałe dynamiki precesji (n_i x dE/dn_i = 0),
(2) w w obu punktach wobec P_S stanów T0' i S' z diagonalizacji QM.
"""
import numpy as np

K = 1.0


def energy(n1, n2, Z):
    return K * n1 @ n2 + Z * np.linalg.norm(n1 + n2) / 2 * (n1[2] - n2[2])


def grad(n1, n2, Z, h=1e-7):
    g = []
    for which in (0, 1):
        gi = np.zeros(3)
        for k in range(3):
            d = np.zeros(3)
            d[k] = h
            a = (n1 + d, n2) if which == 0 else (n1, n2 + d)
            b = (n1 - d, n2) if which == 0 else (n1, n2 - d)
            gi[k] = (energy(*a, Z) - energy(*b, Z)) / (2 * h)
        g.append(gi)
    return g


def config(beta):
    # symetryczne odchylenie od osi x (oś T0 prostopadła do B) o ±beta ku ±z
    n1 = np.array([np.cos(beta), 0.0, np.sin(beta)])
    n2 = np.array([np.cos(beta), 0.0, -np.sin(beta)])
    return n1, n2


print("   x    |  T0': w klas.   P_S QM   |  S': w klas.   P_S QM   | max|n x dE/dn|")
for x in (0.05, 0.1, 0.2756, 0.5, 0.8268, 1.5, 3.0, 10.0):
    Z = x * K
    two_beta = np.arctan(x)
    vals, vecs = np.linalg.eigh(np.array([[-3 * K, 2 * Z], [2 * Z, K]]))
    ps_qm = {"S'": vecs[0, 0] ** 2, "T0'": vecs[0, 1] ** 2}
    row, res = [], 0.0
    for name, b2 in (("T0'", two_beta), ("S'", two_beta - np.pi)):
        n1, n2 = config(b2 / 2)
        g1, g2 = grad(n1, n2, Z)
        res = max(res, np.linalg.norm(np.cross(n1, g1)), np.linalg.norm(np.cross(n2, g2)))
        w = (1 - n1 @ n2) / 2
        row.append("%11.6f  %9.6f" % (w, ps_qm[name]))
    print("%6.4f  |%s   |%s   | %.1e" % (x, row[0], row[1], res))
print("\nx = 0,2756 to Ps 1S przy 1 T, x = 0,8268 przy 3 T (x = 2 g mu_B B / dE_HFS).")
print("Czysty singlet zamiast S' zawyża tempo rozpadu o 1/P_S(S') - 1:")
for x, B in ((0.2756, 1), (0.8268, 3)):
    vals, vecs = np.linalg.eigh(np.array([[-3, 2 * x], [2 * x, 1.0]]))
    print("  %d T: P_S(S') = %.4f  ->  +%.1f %%" % (B, vecs[0, 0] ** 2, 100 * (1 / vecs[0, 0] ** 2 - 1)))

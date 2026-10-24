"""Audit 628: QM first-order alpha^4 reference for the spin-free Hamiltonian of the engine (relativistic kinetic + Coulomb + Darwin
orbit-orbit), against the EBK levels printed by p628.cpp (MODE=levels).

QM with exact hydrogenic expectation values (a = hbar^2/(mu k), E0 = -mu k^2/(2 hbar^2 n^2)):
  <p^4> = 4 mu^2 (E0^2 + 2 E0 k <1/r> + k^2 <1/r^2>)                     (p^2 psi = 2 mu (E0 + k/r) psi)
  kinetic: -(1/m1^3 + 1/m2^3) <p^4>/(8 c^2)
  orbit-orbit, Hermitian ordering p_i (delta_ij + n_i n_j)/r p_j:  H_D = (k q1q2/(2 c^2 m1 m2 |q1q2|/e^2 ...)) -> written with K = k_engine
     <...> = hbar^2 int (|grad psi|^2 + |d_r psi|^2)/r = 2 Re<p^2/r> - hbar^2 l(l+1) <1/r^3> - 4 pi hbar^2 |psi(0)|^2 (last term l = 0 only),
     Re<p^2/r> = 2 mu (E0 <1/r> + k <1/r^2>)
  Darwin contact of each Dirac particle: (pi hbar^2 k/(2 c^2)) (1/m1^2 + 1/m2^2) |psi(0)|^2 (l = 0 only)
  Ferrell (Ps): m alpha^4 c^2 [11/(64 n^4) - 1/(2 n^3 (2l + 1))]  (spin-independent part).
Variants for l = 0: QM_A without any contact term, QM_B with the ordering contact (-4 pi |psi(0)|^2), QM_C = QM_B + Darwin contact.
usage: python3 p628qm.py levels_file [levels_file ...]
"""
import math
import sys


def load(path):
    const, levels = {}, []
    for line in open(path):
        w = line.split()
        if w and w[0] == "CONST":
            const = {w[i]: float(w[i + 1]) for i in range(1, len(w), 2)}
        elif w and w[0] == "LEVEL":
            levels.append((int(w[2]), int(w[4]), float(w[6]), float(w[8]), float(w[10])))
    return const, levels


def qm(const, n, l):
    m1, m2, mu, k, hb, c = const["m1"], const["m2"], const["mu"], const["k"], const["hbar"], const["c"]
    a = hb * hb / (mu * k)
    E0 = -mu * k * k / (2 * hb * hb * n * n)
    r1 = 1 / (n * n * a)
    r2 = 1 / (n ** 3 * (l + 0.5) * a * a)
    r3 = 1 / (n ** 3 * l * (l + 0.5) * (l + 1) * a ** 3) if l > 0 else 0.0
    psi0 = 1 / (math.pi * n ** 3 * a ** 3) if l == 0 else 0.0
    p4 = 4 * mu * mu * (E0 * E0 + 2 * E0 * k * r1 + k * k * r2)
    kin = -(1 / m1 ** 3 + 1 / m2 ** 3) * p4 / (8 * c * c)
    rep2 = 2 * mu * (E0 * r1 + k * r2)
    oo = 2 * rep2 - hb * hb * l * (l + 1) * r3
    pref = -k / (2 * c * c * m1 * m2)          # k q1 q2 = -k for opposite unit charges
    dA = kin + pref * oo
    dB = dA + pref * (-4 * math.pi * hb * hb * psi0)
    dC = dB + math.pi * hb * hb * k / (2 * c * c) * (1 / m1 ** 2 + 1 / m2 ** 2) * psi0
    alpha = k / (hb * c)
    ferrell = m1 * alpha ** 4 * c * c * (11 / (64 * n ** 4) - 1 / (2 * n ** 3 * (2 * l + 1))) if abs(m1 - m2) < 1e-9 * m1 else None
    return dA, dB, dC, ferrell, alpha


for path in sys.argv[1:]:
    const, levels = load(path)
    print("== %s (m1 = %.6e, m2 = %.6e)" % (path, const["m1"], const["m2"]))
    worst = 0.0
    for n, l, E0, E, dE in levels:
        dA, dB, dC, fer, alpha = qm(const, n, l)
        row = "  n=%d l=%d dE_EBK %+.9e J | EBK/QM_A %.7f" % (n, l, dE, dE / dA)
        if l == 0:
            row += "  EBK/QM_B %.7f  EBK/QM_C %.7f" % (dE / dB, dE / dC)
        else:
            worst = max(worst, abs(dE / dA - 1))
        if fer is not None:
            row += "  EBK/Ferrell %.7f  QM_C/Ferrell %.7f" % (dE / fer, dC / fer)
        print(row)
    print("  H628a (l >= 1): max |dE_EBK/dE_QM - 1| = %.3e (alpha^2 = %.3e): %s" % (worst, alpha ** 2, "TAK" if worst < 1e-3 else "NIE"))

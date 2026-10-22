"""Tempa E1 wodoru: reguła punktu środkowego (reguła B) wobec ścisłego QM.

Bez silnika, tylko numpy. Jednostki atomowe, jądro nieskończenie ciężkie.

Reguły klasyczne (foton = jeden mod kołowy klasycznego dipola, import 2):
  * mod współbieżny c+ dla Δl = -1, przeciwbieżny c- dla Δl = +1,
  * tempo A = (2/3) |c_s|^2 ω^3 / c^3   (moc modu / ħω), s = n - n'.
Warianty (orbita, na której liczy się składową Fouriera c_s; częstość ω):
  R0  orbita początkowa n, L = l + 1/2 (Langer); ω = s/n^3 (klasyczna)
  RB  punkt środkowy n_m = (n + n')/2, L_m = max(l, l') (reguła B); ω = s/n_m^3
  RBe punkt środkowy jak RB; ω = dokładna częstość przejścia
  RBw RB razy waga kanału 2 l_>/(2l+1) = L_m/L_i (Langer); dla s -> p waga 2:
      jedyny dozwolony mod (przeciwbieżny) bierze cały udział, bo mod
      współbieżny nie mieści się w siatce (l' = -1)
Dla stanów s porównanie także z sumą modów przeciwbieżnych po wszystkich
harmonicznych na orbicie n o stałym L (jak tempo silnika dla l >= 1, ale z
modem przeciwbieżnym).
QM: A = (4/3) ω^3/c^3 · l_>/(2l+1) · R^2, R — całka radialna.
"""
import numpy as np

C = 137.035999
AU_RATE = 4.134137333e16  # 1/s


def fourier_components(n_orb, L, s, N=8192):
    """c_s^+ i c_s^- dipola x + iy orbity Keplera (a = n^2, L w jedn. ħ)."""
    a = n_orb ** 2
    e = np.sqrt(max(0.0, 1.0 - (L / n_orb) ** 2))
    b = a * np.sqrt(1.0 - e * e)
    u = (np.arange(N) + 0.5) * 2 * np.pi / N
    M = u - e * np.sin(u)
    z = a * (np.cos(u) - e) + 1j * b * np.sin(u)
    w = (1 - e * np.cos(u)) / N  # dM / 2π
    cp = np.sum(z * np.exp(-1j * s * M) * w)
    cm = np.sum(z * np.exp(+1j * s * M) * w)
    return abs(cp), abs(cm)


def laguerre(k, alpha, x):
    if k == 0:
        return np.ones_like(x)
    l0, l1 = np.ones_like(x), 1 + alpha - x
    for j in range(1, k):
        l0, l1 = l1, ((2 * j + 1 + alpha - x) * l1 - (j + alpha) * l0) / (j + 1)
    return l1


R_GRID = np.linspace(0, 900, 360001)


def radial(n, l):
    r = R_GRID
    f = r ** l * np.exp(-r / n) * laguerre(n - l - 1, 2 * l + 1, 2 * r / n)
    f /= np.sqrt(np.trapezoid(f * f * r * r, r))
    return f


_cache = {}


def R(n, l, n2, l2):
    for key in ((n, l), (n2, l2)):
        if key not in _cache:
            _cache[key] = radial(*key)
    return np.trapezoid(_cache[(n, l)] * _cache[(n2, l2)] * R_GRID ** 3, R_GRID)


def omega(n, n2):
    return 0.5 * (1.0 / n2 ** 2 - 1.0 / n ** 2)


def rate_qm(n, l, n2, l2):
    lg = max(l, l2)
    return 4 / 3 * omega(n, n2) ** 3 / C ** 3 * lg / (2 * l + 1) * R(n, l, n2, l2) ** 2


def rate_rule(n, l, n2, l2, variant):
    s = n - n2
    if variant == "R0":
        n_orb, L = n, l + 0.5
    else:
        n_orb, L = 0.5 * (n + n2), max(l, l2)
    w = omega(n, n2) if variant == "RBe" else s / n_orb ** 3
    cp, cm = fourier_components(n_orb, L, s)
    c = cp if l2 == l - 1 else cm
    weight = 2 * max(l, l2) / (2 * l + 1) if variant == "RBw" else 1.0
    return weight * 2 / 3 * c ** 2 * w ** 3 / C ** 3


def counter_mode_sum(n, L, smax=60):
    tot = 0.0
    for s in range(1, smax + 1):
        c = fourier_components(n, L, s, N=16384)[1]
        tot += 2 / 3 * c ** 2 * (s / n ** 3) ** 3 / C ** 3
    return tot


def channels(n, l):
    for n2 in range(1, n):
        for l2 in (l - 1, l + 1):
            if 0 <= l2 < n2:
                yield n2, l2


def lifetime_ns(n, l, fn):
    tot = sum(fn(n, l, n2, l2) for n2, l2 in channels(n, l))
    return 1e9 / (tot * AU_RATE)


if __name__ == "__main__":
    print("Kontrola: A(2p->1s) QM = %.4e 1/s (tablice 6.265e8)"
          % (rate_qm(2, 1, 1, 0) * AU_RATE))
    names = "spdfg"
    V = ("R0", "RB", "RBe", "RBw")
    print("\nCzasy życia [ns]; w nawiasach reguła/QM")
    print("stan        QM   " + "".join("%18s" % v for v in V))
    for l in range(0, 4):
        for n in range(max(l + 2, 3) if l == 0 else l + 2, 9):
            q = lifetime_ns(n, l, rate_qm)
            row = "%d%s  %10.3f  " % (n, names[l], q)
            for v in V:
                t = lifetime_ns(n, l, lambda *a, v=v: rate_rule(*a, v))
                row += "%10.3f (%.3f)" % (t, t / q)
            print(row)
    print("\nStany s: czas życia reguła/QM")
    print("stan   suma modów przeciwbieżnych L=1/2   L=1     RBw")
    for n in range(3, 9):
        q = lifetime_ns(n, 0, rate_qm)
        a = 1e9 / (counter_mode_sum(n, 0.5) * AU_RATE)
        b = 1e9 / (counter_mode_sum(n, 1.0) * AU_RATE)
        rb = lifetime_ns(n, 0, lambda *x: rate_rule(*x, "RBw"))
        print("%ds     %25.3f   %6.3f  %6.3f" % (n, a / q, b / q, rb / q))
    print("\nPojedyncze kanały, tempo reguła/QM:")
    print("kanał      " + "".join("%8s" % v for v in V))
    for (n, l, n2, l2) in [(2, 1, 1, 0), (3, 0, 2, 1), (3, 1, 1, 0), (3, 1, 2, 0),
                           (3, 2, 2, 1), (4, 0, 2, 1), (4, 0, 3, 1), (4, 1, 3, 2),
                           (5, 0, 4, 1), (5, 0, 2, 1), (6, 0, 5, 1), (6, 0, 2, 1),
                           (8, 0, 7, 1), (8, 0, 2, 1)]:
        q = rate_qm(n, l, n2, l2)
        print("%d%s -> %d%s " % (n, names[l], n2, names[l2])
              + "".join("%8.3f" % (rate_rule(n, l, n2, l2, v) / q) for v in V))

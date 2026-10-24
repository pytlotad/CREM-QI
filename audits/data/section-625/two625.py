"""Audit 625: the E1 "against" branches of the spin-1 triplet (611) as two steps through an intermediate state (change of L, then
recoupling with the spin), each at the rule-B midpoint -- variants V1-V5 of prereg625.txt against the exact QM share from the 6j symbol.
Exact QM in rational arithmetic (Racah formula, squared 6j as a Fraction)."""
import math
from fractions import Fraction as F


def fac(n):
    return math.factorial(n)


def tri(a, b, c):
    return abs(a - b) <= c <= a + b


def d2(a, b, c):
    return F(fac(a + b - c) * fac(a - b + c) * fac(-a + b + c), fac(a + b + c + 1))


def sixj2(a, b, c, d, e, f):
    if not (tri(a, b, c) and tri(a, e, f) and tri(d, b, f) and tri(d, e, c)):
        return F(0)
    pre = d2(a, b, c) * d2(a, e, f) * d2(d, b, f) * d2(d, e, c)
    s = 0
    for t in range(max(a + b + c, a + e + f, d + b + f, d + e + c), min(a + b + d + e, b + c + e + f, c + a + f + d) + 1):
        s += (-1) ** t * F(fac(t + 1), fac(t - a - b - c) * fac(t - a - e - f) * fac(t - d - b - f) * fac(t - d - e - c)
                           * fac(a + b + d + e - t) * fac(b + c + e + f - t) * fac(c + a + f + d - t))
    return pre * s * s


def p_qm(l, J, lp, Jp, S=1):
    w = {jp: (2 * jp + 1) * sixj2(l, J, S, jp, lp, 1) for jp in range(abs(lp - S), lp + S + 1)}
    return w[Jp] / sum(w.values())


def cos_tri(a, b, c):
    """cosine of the angle between sides a and b of the triangle (a, b, c)."""
    return (a * a + b * b - c * c) / (2 * a * b)


def born(x):
    return math.sqrt(x * (x + 1))


def langer(x):
    return x + 0.5


def variants(l, J, lp, Jp):
    out = {}
    # V1 / V2: one rotation of the L axis, triangle (L, L', k)
    for name, length, k in (("V1", langer, 1.5), ("V2", born, math.sqrt(2))):
        c = cos_tri(length(l), length(lp), k)
        out[name] = ((1 - c) / 2) ** 2
    # V3 / V4: two steps through mu'' = 0, angles from (l, l', 1) and (J, J', 1)
    for name, length, k in (("V3", born, math.sqrt(2)), ("V4", langer, 1.5)):
        cl = cos_tri(length(l), length(lp), k)
        cj = cos_tri(length(J), length(Jp), k)
        out[name] = (1 - cl * cl) * (1 - cj * cj) / 4
    # V5: both rotations composed about one axis (all intermediate states)
    cl = cos_tri(born(l), born(lp), math.sqrt(2))
    cj = cos_tri(born(J), born(Jp), math.sqrt(2))
    th = math.acos(cl) + math.acos(cj)
    out["V5"] = ((1 - math.cos(th)) / 2) ** 2
    return out


families = (("F1: l -> l-1, J = l-1 -> J' = l", lambda l: (l, l - 1, l - 1, l)),
            ("F2: l -> l+1, J = l+1 -> J' = l", lambda l: (l, l + 1, l + 1, l)))
LS = (2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 100, 300, 1000)
names = ("V1", "V2", "V3", "V4", "V5")
summary = {}
for fam, args in families:
    print(fam)
    print("   l      p_QM          l^4 p_QM   " + "   ".join("%s/QM" % n for n in names) + "   l^5 (p_V3 - p_QM)")
    for l in LS:
        a = args(l)
        q = p_qm(*a)
        v = variants(*a)
        r = {n: v[n] / float(q) for n in names}
        for n in names:
            summary.setdefault((fam, n), []).append((l, r[n]))
        print("  %4d  %.6e  %.6f  " % (l, float(q), float(q) * l ** 4)
              + "  ".join("%.6f" % r[n] for n in names)
              + "   %+.4f" % ((v["V3"] - float(q)) * l ** 5))
    # exact rational form for small l (pattern)
    print("   p_QM dokladnie: " + ", ".join("l=%d: %s" % (l, p_qm(*args(l))) for l in (2, 3, 4, 5, 6)))
print()
print("REGULY (|p_V/p_QM - 1|):")
for n in names:
    asym = all(abs(dict(summary[(f, n)])[1000] - 1) < 0.01
               and abs(dict(summary[(f, n)])[1000] - 1) < abs(dict(summary[(f, n)])[100] - 1) for f, _ in families)
    exact = max(abs(r - 1) for f, _ in families for l, r in summary[(f, n)] if l <= 40)
    print("  %s: H625a %s (l = 1000: %s), H625b %s (max l <= 40: %.3e)" % (
        n, "TAK" if asym else "NIE", " / ".join("%.5f" % (dict(summary[(f, n)])[1000] - 1) for f, _ in families),
        "TAK" if exact < 1e-9 else "NIE", exact))

# POST HOC (after the result above, not a pre-registered variant): the exact shares have the closed forms
# F1 p = 1/(l (2l - 1))^2, F2 p = 1/((l + 1)(2l + 3))^2, i.e. p = [ (1/2) (1/l_>) (1/(J + 1/2)) ]^2 with J the INITIAL J:
# step 1 at the rule-B midpoint (sin theta_L = hbar/L_m, L_m = l_>, exact for Born lengths), step 2 with the rule-Q length of the initial J.
print()
print("POST HOC: p_QM wobec [(1/2)(1/l_>)(1/(J + 1/2))]^2 i wobec postaci zamknietych, l = 2 ... 200 (dokladnie, Fraction):")
for fam, args in families:
    worst_closed = worst_form = 0
    for l in range(2, 201):
        lo, J, lp, Jp = args(l)
        q = p_qm(lo, J, lp, Jp)
        lg = max(lo, lp)
        closed = F(1, (l * (2 * l - 1)) ** 2) if lp == l - 1 else F(1, ((l + 1) * (2 * l + 3)) ** 2)
        form = (F(1, 2) * F(1, lg) * F(2, 2 * J + 1)) ** 2
        worst_closed = max(worst_closed, abs(q - closed))
        worst_form = max(worst_form, abs(q - form))
    print("  %s: max |p_QM - postac zamknieta| = %s, max |p_QM - forma| = %s" % (fam, worst_closed, worst_form))

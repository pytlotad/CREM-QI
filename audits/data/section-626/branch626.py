"""Audit 626: the length pair of 625d (orbital step at the rule-B midpoint, spin step at the rule-Q length of the initial J) on all E1
branches of the 6j symbol -- along, sideways, against -- for spin 1/2 (hydrogen) and spin 1 (Ps triplet); variants W1-W3 of prereg626.txt.
Exact QM from the Racah formula in rational arithmetic (half-integers carried as doubled integers)."""
import math
from fractions import Fraction as F


def fac2(n2):
    """factorial of n2/2 for even n2 (all 6j arguments combine to integers)."""
    assert n2 % 2 == 0, n2
    return math.factorial(n2 // 2)


def tri2(a, b, c):
    return abs(a - b) <= c <= a + b and (a + b + c) % 2 == 0


def d2(a, b, c):
    return F(fac2(a + b - c) * fac2(a - b + c) * fac2(-a + b + c), fac2(a + b + c + 2))


def sixj2(a, b, c, d, e, f):
    """squared 6j with all arguments doubled."""
    if not (tri2(a, b, c) and tri2(a, e, f) and tri2(d, b, f) and tri2(d, e, c)):
        return F(0)
    pre = d2(a, b, c) * d2(a, e, f) * d2(d, b, f) * d2(d, e, c)
    s = 0
    lo = max(a + b + c, a + e + f, d + b + f, d + e + c)
    hi = min(a + b + d + e, b + c + e + f, c + a + f + d)
    for t2 in range(lo, hi + 1, 2):
        s += (-1) ** (t2 // 2) * F(fac2(t2 + 2), fac2(t2 - a - b - c) * fac2(t2 - a - e - f) * fac2(t2 - d - b - f)
                                    * fac2(t2 - d - e - c) * fac2(a + b + d + e - t2) * fac2(b + c + e + f - t2)
                                    * fac2(c + a + f + d - t2))
    return pre * s * s


def qm_branches(l, J2, lp, S2):
    """QM shares J -> J' (J', S doubled) for l -> lp."""
    w = {}
    for Jp2 in range(abs(2 * lp - S2), 2 * lp + S2 + 1, 2):
        w[Jp2] = (Jp2 + 1) * sixj2(2 * l, J2, S2, Jp2, 2 * lp, 2)
    z = sum(w.values())
    return {k: v / z for k, v in w.items() if z > 0}


def dsq(S2, m2p, m2, th):
    """|d^S_{m' m}(theta)|^2 for S = 1/2 and S = 1 (projections doubled)."""
    c, s = math.cos(th), math.sin(th)
    if S2 == 1:
        return math.cos(th / 2) ** 2 if m2p == m2 else math.sin(th / 2) ** 2
    if abs(m2p) > 2 or abs(m2) > 2:
        return 0.0
    table = {(2, 2): ((1 + c) / 2) ** 2, (-2, -2): ((1 + c) / 2) ** 2, (2, -2): ((1 - c) / 2) ** 2, (-2, 2): ((1 - c) / 2) ** 2,
             (0, 0): c * c}
    return table.get((m2p, m2), s * s / 2)


def model(l, J2, lp, S2, variant):
    dl = lp - l
    mu2 = J2 - 2 * l
    sinL = 1.0 / max(l, lp)                     # rule-B midpoint, Born lengths: sin theta_L = 1/l_>
    thL = math.asin(sinL)
    out = {}
    jps = list(range(abs(2 * lp - S2), 2 * lp + S2 + 1, 2))
    for Jp2 in jps:
        dJ2 = Jp2 - J2
        if abs(dJ2) > 2:
            out[Jp2] = 0.0
            continue
        mup2 = Jp2 - 2 * lp
        if variant == "W3":
            sinJ = 1.0 / (max(J2, Jp2) / 2.0) if dJ2 != 0 else 0.0
        else:
            sinJ = 2.0 / (J2 + 1)                 # 1/(J + 1/2) of the initial J
        thJ = math.asin(min(1.0, sinJ))
        if variant == "W1":
            cosT = 1.0 - sinL * (2.0 / (J2 + 1))
            out[Jp2] = dsq(S2, mup2, mu2, math.acos(max(-1.0, cosT)))
            continue
        kind = "along" if dJ2 == 2 * dl else ("side" if dJ2 == 0 else "against")
        if kind == "side":
            out[Jp2] = dsq(S2, mup2, mu2, thL)
        elif kind == "against":
            mid2 = mu2 - 2 * dl
            out[Jp2] = dsq(S2, mid2, mu2, thL) * dsq(S2, mup2, mid2, thJ) if abs(mid2) <= S2 else 0.0
        else:
            out[Jp2] = None
    if variant != "W1":
        rest = 1.0 - sum(v for v in out.values() if v is not None)
        for k, v in out.items():
            if v is None:
                out[k] = rest
    return out


def kind_of(J2, Jp2, dl):
    dJ2 = Jp2 - J2
    return "along" if dJ2 == 2 * dl else ("side" if dJ2 == 0 else ("against" if dJ2 == -2 * dl else "forbidden"))


variants = ("W1", "W2", "W3")
worst = {(v, S2, k): (0.0, None) for v in variants for S2 in (1, 2) for k in ("along", "side", "against")}
scaling = {}
LMAX = 20
for S2, tag in ((1, "wodor, S = 1/2"), (2, "tryplet Ps, S = 1")):
    print("== %s" % tag)
    for l in list(range(0, LMAX + 1)) + [40, 80]:
        for dl in (-1, 1):
            lp = l + dl
            if lp < 0:
                continue
            for J2 in range(abs(2 * l - S2), 2 * l + S2 + 1, 2):
                qm = qm_branches(l, J2, lp, S2)
                if not qm:
                    continue
                for v in variants:
                    m = model(l, J2, lp, S2, v)
                    for Jp2, q in qm.items():
                        k = kind_of(J2, Jp2, dl)
                        if k == "forbidden":
                            continue
                        dev = abs(m[Jp2] - float(q))
                        if l in (10, 20, 40, 80):
                            key = (v, S2, k, l)
                            scaling[key] = max(scaling.get(key, 0.0), dev)
                        if l <= LMAX and dev > worst[(v, S2, k)][0]:
                            worst[(v, S2, k)] = (dev, (l, lp, J2 / 2, Jp2 / 2, m[Jp2], float(q)))
                if l in (1, 2, 3) and S2 == 2 or l in (0, 1, 2) and S2 == 1:
                    m2 = model(l, J2, lp, S2, "W2")
                    print("  l=%d -> %d, J=%g: " % (l, lp, J2 / 2) + "  ".join(
                        "J'=%g [%s] QM %.6f W2 %.6f" % (Jp2 / 2, kind_of(J2, Jp2, dl), float(q), m2[Jp2]) for Jp2, q in sorted(qm.items())))
print()
print("H626a / H626b: max |P_V - P_QM| (l <= %d) i gdzie; skalowanie przy l = 10 / 20 / 40 / 80" % LMAX)
for v in variants:
    for S2 in (1, 2):
        for k in ("along", "side", "against"):
            dev, where = worst[(v, S2, k)]
            sc = " / ".join("%.2e" % scaling.get((v, S2, k, l), 0.0) for l in (10, 20, 40, 80))
            print("  %s S=%s %-7s max %.3e %s | l=10/20/40/80: %s" % (
                v, "1/2" if S2 == 1 else "1", k, dev,
                "(l=%d->%d J=%g J'=%g: %.6f / %.6f)" % where if where else "", sc))
    exact = max(worst[(v, S2, k)][0] for S2 in (1, 2) for k in ("along", "side", "against"))
    print("  %s: H626a %s (max %.3e)" % (v, "TAK" if exact < 1e-12 else "NIE", exact))

# POST HOC (after the result above; not a pre-registered variant).  The sideways shares have closed forms (exact, l = 2 ... 60):
#   S = 1/2: 1/(l_> (2 l' + 1));  S = 1, mu = 0: 1/l_>^2 = sin^2 theta_L;  S = 1, mu = -+1: (2l + 1)/(l_>^2 (2 l' + 1)).
# Asymptotically each is twice the one-step |d^S(theta_L)|^2 of W2.  W2x: sideways = 2 |d^S_{mu' mu}(theta_L)|^2, along by unitarity.
print()
print("POST HOC 1: postacie zamkniete galezi w bok (dokladnie, l = 2 ... 60):")
for S2 in (1, 2):
    for dl in (-1, 1):
        for J2off in range(-S2, S2 + 1, 2):
            bad, n = 0, 0
            for l in range(2, 61):
                lp, J2 = l + dl, 2 * l + J2off
                qm = qm_branches(l, J2, lp, S2)
                if J2 not in qm:
                    continue
                lg = max(l, lp)
                if S2 == 1:
                    closed = F(1, lg * (2 * lp + 1))
                elif J2off == 0:
                    closed = F(1, lg * lg)
                else:
                    closed = F(2 * l + 1, lg * lg * (2 * lp + 1))
                n += 1
                bad += qm[J2] != closed
            if n:
                print("  S=%s dl=%+d J=l%+g: %d / %d przypadkow rownych postaci zamknietej" % ("1/2" if S2 == 1 else "1", dl, J2off / 2, n - bad, n))


def model_x(l, J2, lp, S2):
    m = model(l, J2, lp, S2, "W2")
    dl = lp - l
    mu2 = J2 - 2 * l
    thL = math.asin(1.0 / max(l, lp))
    out = {}
    for Jp2 in m:
        k = kind_of(J2, Jp2, dl)
        out[Jp2] = 2 * dsq(S2, Jp2 - 2 * lp, mu2, thL) if k == "side" else (m[Jp2] if k == "against" else None)
    rest = 1.0 - sum(v for v in out.values() if v is not None)
    return {k: (rest if v is None else v) for k, v in out.items()}


print("POST HOC 2: W2x (w bok = 2 |d(theta_L)|^2, wbrew = 625d, wzdluz z unitarnosci), max |P - P_QM| bez celu l' = 0, oraz l^k skalowanie:")
for S2 in (1, 2):
    for k in ("along", "side", "against"):
        res = []
        for l in (2, 5, 10, 20, 40, 80, 160):
            dev = 0.0
            for dl in (-1, 1):
                lp = l + dl
                for J2 in range(abs(2 * l - S2), 2 * l + S2 + 1, 2):
                    qm = qm_branches(l, J2, lp, S2)
                    mx = model_x(l, J2, lp, S2)
                    for Jp2, q in qm.items():
                        if kind_of(J2, Jp2, dl) == k:
                            dev = max(dev, abs(mx[Jp2] - float(q)))
            res.append((l, dev))
        print("  S=%s %-7s " % ("1/2" if S2 == 1 else "1", k) + "  ".join("l=%d %.2e (l^3 x %.3f)" % (l, d, d * l ** 3) for l, d in res))

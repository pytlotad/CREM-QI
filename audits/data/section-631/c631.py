"""Audit 631, H631c: o-Ps 3d-start runs with CREM_DEBUG_BORN -- the photon draw (J -> J'), label after the turn, J' histogram against 6j,
and the next photon (2p -> 1s) energy grouped by J'.  usage: python3 c631.py file.err [file.err ...]"""
import math
import re
import sys
from fractions import Fraction as F

EV = 1.602176634e-19


def fac(n): return math.factorial(n)
def tri(a, b, c): return abs(a - b) <= c <= a + b
def d2(a, b, c): return F(fac(a + b - c) * fac(a - b + c) * fac(-a + b + c), fac(a + b + c + 1))


def sixj2(a, b, c, d, e, f):
    if not (tri(a, b, c) and tri(a, e, f) and tri(d, b, f) and tri(d, e, c)):
        return F(0)
    pre = d2(a, b, c) * d2(a, e, f) * d2(d, b, f) * d2(d, e, c)
    s = 0
    for t in range(max(a + b + c, a + e + f, d + b + f, d + e + c), min(a + b + d + e, b + c + e + f, c + a + f + d) + 1):
        s += (-1) ** t * F(fac(t + 1), fac(t - a - b - c) * fac(t - a - e - f) * fac(t - d - b - f) * fac(t - d - e - c)
                            * fac(a + b + d + e - t) * fac(b + c + e + f - t) * fac(c + a + f + d - t))
    return pre * s * s


runs = []
for path in sys.argv[1:]:
    cur = None
    for line in open(path):
        if line.startswith("START idx"):
            cur = {"idx": int(line.split()[2]), "photon": None, "eg": []}
            runs.append(cur)
        elif cur is None:
            continue
        elif line.startswith("BORN photon"):
            m = re.search(r"l=(\d+) J=(\d+) -> l'=(\d+) J'=(\d+) \(basin (\d+), after turn (\d+)\)", line)
            cur["photon"] = tuple(int(x) for x in m.groups())
        elif line.startswith("CREM_LEDGER t"):
            cur["eg"].append(float(line.split()[4]) / EV)
withdraw = [r for r in runs if r["photon"]]
print("przebiegow %d, z losowaniem przy fotonie %d" % (len(runs), len(withdraw)))
bad = [r["idx"] for r in withdraw if r["photon"][3] != r["photon"][5]]
print("etykieta po obrocie != J': %d %s" % (len(bad), bad))
changed = sum(1 for r in withdraw if r["photon"][3] != r["photon"][4])
print("J' inne niz basen klasyczny: %d z %d" % (changed, len(withdraw)))
# histogram (J, J') and chi^2 against 6j
hist = {}
for r in withdraw:
    l, J, lp, Jp = r["photon"][:4]
    hist.setdefault((l, J, lp), {}).setdefault(Jp, 0)
    hist[(l, J, lp)][Jp] += 1
chi, dof = 0.0, 0
for (l, J, lp), h in sorted(hist.items()):
    n = sum(h.values())
    w = {Jp: (2 * Jp + 1) * sixj2(l, J, 1, Jp, lp, 1) for Jp in range(abs(lp - 1), lp + 2)}
    z = sum(w.values())
    exp = {k: float(v / z) * n for k, v in w.items()}
    print("  l=%d J=%d -> l'=%d: n=%d  J': %s  (6j: %s)" % (l, J, lp, n, " ".join("%d:%d" % (k, h.get(k, 0)) for k in sorted(exp)),
                                                     " ".join("%d:%.2f" % kv for kv in sorted(exp.items()))))
    merged_o = merged_e = 0
    for k in sorted(exp, key=lambda k: exp[k]):
        merged_o += h.get(k, 0)
        merged_e += exp[k]
        if merged_e >= 5:
            chi += (merged_o - merged_e) ** 2 / merged_e
            dof += 1
            merged_o = merged_e = 0
    dof -= 1
print("chi^2 histogramu = %.2f przy %d stopniach swobody (raport)" % (chi, max(dof, 0)))
# next photon energy by J'
groups = {}
for r in withdraw:
    if len(r["eg"]) >= 2:
        groups.setdefault(r["photon"][3], []).append(r["eg"][1])
for Jp, es in sorted(groups.items()):
    print("  J'=%d: foton 2p -> 1s %d x, srednia %.9f eV, rozrzut %.2e eV" % (Jp, len(es), sum(es) / len(es), max(es) - min(es)))

# Audit 782 POST HOC (does not change verdicts): (1) total Fermi rate out of each start mode (header #F of P10 file) vs threshold 2 Gamma (weight damping of an undriven mode
# relative to the driven one); (2) log-slope of lower-mode weights (#W) per period in units of GT; (3) radiated energy vs Fermi expectation.
import os, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); GT = 3e-4
pr = lambda *a: print(*a, flush=True)
LN = [(0,0),(0,1),(0,-1),(1,0),(0,2),(0,-2),(1,1),(1,-1),(2,0),(1,2),(1,-2),(2,1),(2,-1)]
tot, F = {}, {}
for l in open(os.path.join(O, 'P10_0_1.txt')):
    if l.startswith('#F'):
        p = l.split(); a = (int(p[1]), int(p[2])); b = (int(p[4]), int(p[5])); g = float(p[11]); tot[a] = tot.get(a, 0) + g; F[(a, b)] = g
pr("(1) suma temp Fermiego z modu / Gamma przy DIP 10 (DIP 1: /10); prog wzrostu nizszego modu ~2 (tlumienie wagi modu nienapedzanego 2 Gamma):")
for a in [(0, 1), (0, 2), (1, 0), (1, 1), (1, -2), (2, 0)]:
    pr(f"  {str(a):8s} DIP 10: {tot[a]:7.3f}   DIP 1: {tot[a] / 10:6.3f}   najsilniejszy kanal: " + max(((g, b) for (x, b), g in F.items() if x == a), key=lambda z: z[0]).__repr__())
def W(f):
    a = np.array([[float(x) for x in l.split()[1:]] for l in open(os.path.join(O, f)) if l.startswith('#W')]); return a[:, 0], a[:, 1:]
pr("(2) nachylenie ln w(mod nizszy) na okres, w jednostkach GT; przewidywanie liniowe Gamma_F(kanal)/Gamma - 2:")
for f, col, lo, hi, ch in [('P10_0_2.txt', 1, 3000, 9600, ((0, 2), (0, 1))), ('P1_0_2.txt', 1, 3000, 9600, ((0, 2), (0, 1))),
                           ('P10_1_1.txt', 0, 2000, 9600, ((1, 1), (0, 0))), ('P1_1_1.txt', 0, 2000, 9600, ((1, 1), (0, 0))),
                           ('P10_0_1.txt', 0, 400, 1200, ((0, 1), (0, 0))), ('P1_0_1.txt', 0, 400, 1200, ((0, 1), (0, 0)))]:
    k, w = W(f); s = (k >= lo) & (k <= hi); c = np.polyfit(k[s], np.log(w[s, col] + 1e-30), 1); sc = 10 if f.startswith('P10') else 1
    pr(f"  {f:12s} {str(LN[col]):8s} okresy {lo}-{hi}: {c[0] / GT:+7.2f} GT   przewidywanie {F[ch] * sc / 10 - 2:+6.2f}")
pr("(3) energia wypromieniowana na okres (mediana ost. 1000) wobec Gamma_F,tot * (E_start - E_(0,0)) * Tr0 (gdyby stan sie rozpadal):")
B = 0.2; En = lambda t: -0.5 / (t[0] + 0.5 + np.sqrt(t[1] ** 2 + 2 * B)) ** 2; Tr0 = 2 * np.pi * (0.5 + np.sqrt(1.4)) ** 3
for r in (1, 10):
    for a in [(1, 1), (1, -2), (2, 0), (1, 0)]:
        q = np.loadtxt(os.path.join(O, f"P{r}_{a[0]}_{a[1]}.txt"), comments='#', ndmin=2); ex = tot[a] * r / 10 * GT * (En(a) - En((0, 0)))
        pr(f"  DIP {r:2d} {str(a):8s} E_rad {np.median(q[-1000:, 8]):.2e}  oczekiwane przy rozpadzie ~{ex:.2e}  stosunek {np.median(q[-1000:, 8]) / ex:.1e}")

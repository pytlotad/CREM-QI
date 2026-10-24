"""Trzy sprawdzenia do rekomendacji po audycie 614 (bez silnika, numpy).

Skrypt autora (dolaczony 2026-10-24), skopiowany bez zmian poza sciezka modulu
(sys.path -> ../section-560, gdzie lezy tempo_e1_punkt_srodkowy.py).

1. Kontakt s jako granica l -> 0 wzorów reguły B dla l >= 1 (czy #10 jest
   potrzebny dla nadsubtelnej). Rachunek dokładny na ułamkach, g = 2(1 + a).
   Wodór: A_lj ∝ <r^-3>[2<L·J> + g(3(S·r)^2 - <S·J>)]/(j(j+1)), j = l + 1/2,
   <r^-3> ∝ 1/(l(l+1/2)(l+1)), (S·r)^2 = 1/4 (reguła 521); Fermi ∝ (16/3)(g/2).
   Ps: tryplet J = l + 1, a<L·S> + b<S12> (import 18, a ∝ (2g-1)/2, b ∝ (g/2)^2/4),
   singlet l >= 1 bez członów spinowych; kontakt magnetyczny T - S = (g/2)^2 · 8/3.
   Jednostki: wodór (μ0/4π) μ_B g_p μ_N/(a³n³), Ps αħ³/(m²c a³n³).
2. Tempo E1 z iloczynu amplitud na obu końcach (|c|² -> c_i c_f, reguła B 594)
   wobec QM — czy jedna zasada zastępuje hybrydę 585/613.
3. Logarytm Bethego stanu s z widma klasycznej orbity Langera: średnia ln(2sω)
   z wagą |p_s|²ω_s ∝ s³|x_s|² (najprostsza korespondencja, bez różnic reguły B).
"""
import sys
from fractions import Fraction as Fr
import numpy as np

sys.path.insert(0, '../section-560')
import tempo_e1_punkt_srodkowy as T


def hydrogen_continuation(a):
    g = 2 * (1 + a)
    l, j = 0, Fr(1, 2)
    # granica: l<L·J>/(l(l+1/2)(l+1)) -> (2l+3)/(2(l+1/2)(l+1)); spin -l/2 -> -1/2/((l+1/2)(l+1))
    orb = 2 * Fr(2 * l + 3, 2) / ((l + Fr(1, 2)) * (l + 1))
    spin = g * Fr(-1, 2) / ((l + Fr(1, 2)) * (l + 1))
    return (orb + spin) / (j * (j + 1)), Fr(16, 3) * g / 2


def ps_continuation(a):
    g = 2 * (1 + a)
    so = (2 * g - 1) / 2 * 2          # a<L·S>: l/(l(l+1/2)(l+1)) -> 2
    ten = (g / 2) ** 2 / 4 * Fr(-4, 3)  # b<S12>: -2l/((2l+3) l(l+1/2)(l+1)) -> -4/3
    return so + ten, (g / 2) ** 2 * Fr(8, 3)


def check_contact():
    print("1. Kontakt s jako granica l -> 0 (wartość ciągła / kontakt):")
    for name, f in (("wodór 1s, A", hydrogen_continuation), ("Ps 1S, T - S magnet.", ps_continuation)):
        c0, k0 = f(Fr(0))
        print("   %-22s g = 2: %s / %s" % (name, c0, k0))
        a = Fr(115965218, 10 ** 11)
        c, k = f(a)
        print("   %-22s g zmierzone: różnica względna %.2e" % ("", float((c - k) / k)))
    c, k = hydrogen_continuation(Fr(115965218, 10 ** 11))
    print("   wodór 1s: %.2f MHz z 1420,406 MHz" % (1420.406 * float((k - c) / k)))


def rate_endpoint(n, l, n2, l2, weighted):
    s = n - n2
    sense = 0 if l2 == l - 1 else 1
    ci = T.fourier_components(n, l + 0.5, s)[sense]
    cf = T.fourier_components(n2, l2 + 0.5, s)[sense]
    wt = 2 * max(l, l2) / (2 * l + 1) if weighted else 1.0
    return wt * 2 / 3 * ci * cf * T.omega(n, n2) ** 3 / T.C ** 3


def check_endpoint():
    print("\n2. Tempo E1 z iloczynu amplitud na końcach, czas życia reguła/QM (z wagą 2l_>/(2l+1)):")
    names = "spdf"
    row = []
    for (n, l) in [(2, 1), (3, 0), (3, 1), (3, 2), (4, 0), (4, 1), (4, 3), (6, 0), (6, 1), (8, 1), (8, 3)]:
        q = T.lifetime_ns(n, l, T.rate_qm)
        r = T.lifetime_ns(n, l, lambda *x: rate_endpoint(*x, True))
        row.append("%d%s %.2f" % (n, names[l], r / q))
    print("   " + ", ".join(row))


def check_bethe():
    print("\n3. Logarytm Bethego stanów s z widma orbity Langera (L = ħ/2):")
    qm = {1: 2.984, 2: 2.812}
    for n in (1, 2):
        w = 1 / n ** 3
        num = den = 0.0
        for s in range(1, 400):
            cp, cm = T.fourier_components(n, 0.5, s, N=65536)
            W = s ** 3 * (cp ** 2 + cm ** 2)
            num += W * np.log(2 * s * w)
            den += W
        print("   %ds: klasycznie %.3f, QM %.3f" % (n, num / den, qm[n]))


if __name__ == "__main__":
    check_contact()
    check_endpoint()
    check_bethe()

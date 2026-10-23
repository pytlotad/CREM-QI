"""Importy QED 11, 12, 16 w rzędzie wiodącym i z poprawką jednopętlową.

Bez silnika. Model (Model.md, sekcja 9): p-Ps 124,495 ps, o-Ps 138,68 ns,
HFS 1s 204,66 GHz (część magnetyczna z g zmierzonym, czyli z (g/2)^2).
Współczynniki poprawek (literatura QED):
  p-Ps:  Γ = Γ0 [1 - (5 - π²/4) α/π + 2 α² ln(1/α)]      (Harris, Brown 1957; Khriplovich, Yelkhovsky 1990)
  o-Ps:  Γ = Γ0 [1 - 10,286606 α/π - (1/3) α² ln α]       (Caswell, Lepage, Sapirstein 1977; Adkins i in.)
  HFS:   ΔE = α⁴mc² [7/12 - (α/π)(8/9 + ½ ln 2)]          (Karplus, Klein 1952)
"""
# Skrypt autora (propozycja 1.4), zachowany bez zmian ponizej docstringu jako zrodlo.
# UWAGA (audyt 575): czlon o-Ps ma w literaturze postac -(1/3) alpha^2 ln(1/alpha);
# "- a*a*log(a)/3" ponizej daje znak przeciwny.  Wariant poprawny: warianty575.py.
import math

a = 7.2973525693e-3
api = a / math.pi
mc2_hbar = 510998.95 / 6.582119569e-16      # 1/s
mc2_h = mc2_hbar / (2 * math.pi)            # Hz

G_p0 = a ** 5 * mc2_hbar / 2
G_o0 = 2 * (math.pi ** 2 - 9) / (9 * math.pi) * a ** 6 * mc2_hbar
G_p1 = G_p0 * (1 - (5 - math.pi ** 2 / 4) * api + 2 * a * a * math.log(1 / a))
G_o1 = G_o0 * (1 - 10.286606 * api - a * a * math.log(a) / 3)

tp_meas, to_meas = 125.14, 142.04  # ps, ns
print("              LO          + jedna pętla     pomiar")
print("tau p-Ps  %8.3f ps   %8.3f ps      %.2f ps  (model dziś 124,495)"
      % (1e12 / G_p0, 1e12 / G_p1, tp_meas))
print("tau o-Ps  %8.3f ns   %8.3f ns      %.2f ns  (model dziś 138,68)"
      % (1e9 / G_o0, 1e9 / G_o1, to_meas))
print("o/p       %8.1f      %8.1f         %.1f    (model dziś 1114)"
      % (G_p0 / G_o0, G_p1 / G_o1, to_meas * 1e3 / tp_meas))
print("odchyłki z poprawką: p %+.3f %%, o %+.3f %%"
      % (100 * (1e12 / G_p1 / tp_meas - 1), 100 * (1e9 / G_o1 / to_meas - 1)))

E4 = a ** 4 * mc2_h / 1e9                  # GHz
hfs_lo = 7 / 12 * E4
magn = E4 / 3
a_e = 0.00115965218
anomaly_in_model = magn * ((1 + a_e) ** 2 - 1)   # model już ma (g/2)^2 w części magnetycznej
one_loop = -api * (8 / 9 + 0.5 * math.log(2)) * E4
model = 204.66
print("\nHFS 1s: LO %.3f GHz; model %.2f (LO + anomalia %.3f)" % (hfs_lo, model, anomaly_in_model))
print("poprawka jednopętlowa bez części anomalnej: %.3f GHz" % (one_loop - anomaly_in_model))
print("model + poprawka: %.3f GHz wobec 203,394 (Ishida i in. 2014)"
      % (model + one_loop - anomaly_in_model))

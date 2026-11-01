"""Audit 682: is the anomalous part of the s-state contact the l -> 0 limit of rule B's l >= 1 formulas?

Offline, exact (sympy).  Formulas as in ../section-622/sprawdzenia_po_614.py (the author's script) and the engine
(pStateHyperfineEnergy, audit 520; Ps a<L.S> + b<S12>, audit 431), g = 2(1 + a).
  hydrogen  A_lj ~ <r^-3> [2<L.J> + g (3 (S.r^)^2 - <S.J>)] / (j(j+1)),  <r^-3> ~ 1/(l(l+1/2)(l+1)),  (S.r^)^2 = 1/4;
            contact (Fermi) ~ (16/3)(g/2)
  Ps        triplet a<L.S> + b<S12>,  a ~ (2g - 1)/2,  b ~ (g/2)^2/4  (times <r^-3>);  magnetic contact T - S ~ (g/2)^2 8/3
C1 limit in l (622), C2 the same in rule B's Langer lengths (L = l + 1/2), C3 the (2j+1)-weighted mean over the j / J
branches, C4 limit of the geometry at g = 2 times the measured moments (today's engine).
"""
import sympy as sp

l, a, L = sp.symbols('l a L')
g = 2 * (1 + a)
half = sp.Rational(1, 2)
r3 = 1 / (l * (l + half) * (l + 1))
A_E = sp.Rational(115965218, 10**11)  # electron anomaly, measured (CODATA)


def h_branch(j):
    jj = j * (j + 1)
    LJ = (jj + l * (l + 1) - sp.Rational(3, 4)) / 2
    SJ = (jj + sp.Rational(3, 4) - l * (l + 1)) / 2
    return r3 * (2 * LJ + g * (sp.Rational(3, 4) - SJ)) / jj


def ps_branch(J):
    # <L.S> and <S12> of the S = 1 triplet with total J
    LS = (J * (J + 1) - l * (l + 1) - 2) / 2
    S12 = {l + 1: -2 * l / (2 * l + 3), l: sp.Integer(2), l - 1: -2 * (l + 1) / (2 * l - 1)}[J]
    return r3 * ((2 * g - 1) / 2 * LS + (g / 2) ** 2 / 4 * S12)


H_contact = sp.Rational(16, 3) * g / 2
P_contact = (g / 2) ** 2 * sp.Rational(8, 3)


def report(name, expr, contact, out):
    lim = sp.simplify(sp.limit(sp.simplify(expr), l, 0))
    if lim in (sp.oo, -sp.oo, sp.zoo) or lim.has(sp.oo):
        out.append(f"  {name:44s} limit diverges: {lim}")
        return None
    at2 = sp.simplify(lim.subs(a, 0))
    dlog = sp.simplify(sp.diff(sp.log(lim), a).subs(a, 0)) if at2 != 0 else None
    rel = float(((lim - contact) / contact).subs(a, A_E))
    out.append(f"  {name:44s} limit {sp.factor(lim)}  | g=2: {at2} (contact {contact.subs(a, 0)})"
               f"  | dln/da {dlog} (contact {sp.diff(sp.log(contact), a).subs(a, 0)})  | rel. gap at a_e {rel:+.4e}")
    return rel


out = []
out.append("C1 -- limit in l (622):")
relH = report("H  j = l + 1/2", h_branch(l + half), H_contact, out)
relP = report("Ps J = l + 1", ps_branch(l + 1), P_contact, out)
out.append("C2 -- the same in Langer lengths, l = L - 1/2, L -> 1/2:")
for name, e, k in (("H  j = l + 1/2", h_branch(l + half), H_contact), ("Ps J = l + 1", ps_branch(l + 1), P_contact)):
    eL = e.subs(l, L - half)
    lim = sp.simplify(sp.limit(sp.simplify(eL), L, half))
    out.append(f"  {name:44s} limit {sp.factor(lim)}  (= C1: {sp.simplify(lim - sp.limit(sp.simplify(e), l, 0)) == 0})")
out.append("C3 -- branches and (2j+1)-weighted mean:")
report("H  j = l - 1/2 alone", h_branch(l - half), H_contact, out)
wH = ((2 * l + 2) * h_branch(l + half) + 2 * l * h_branch(l - half)) / (2 * (2 * l + 1))
report("H  mean over j (weights (2j+1)/(2(2l+1)))", wH, H_contact, out)
for J, nm in ((l, "Ps J = l alone"), (l - 1, "Ps J = l - 1 alone")):
    report(nm, ps_branch(J), P_contact, out)
wP = sum((2 * J + 1) * ps_branch(J) for J in (l + 1, l, l - 1)) / (3 * (2 * l + 1))
out.append(f"  {'Ps mean over J (weights (2J+1)/(3(2l+1)))':44s} = {sp.simplify(wP)} for every l (trace of L.S and S12)")
out.append("C4 -- geometry at g = 2, moments after: contact = [limit at g = 2] x (g/2) (H), x (g/2)^2 (Ps) = Fermi exactly (by construction).")
out.append("")

# Numbers against measurement.  h in eV s, c in m/s.
h, c = 4.135667696e-15, 299792458.0
alpha, mc2 = 1 / 137.035999084, 510998.95
fH = 1420.405751768e6        # Hellwig et al. 1970
fPs = 203.3942e9             # Ishida et al. 2014
out.append(f"Lines: H 1s HFS {fH/1e6:.9f} MHz = {c/fH*1e3:.4f} mm = {h*fH*1e6:.6f} ueV;"
           f"  Ps 1S HFS {fPs/1e9:.4f} GHz = {c/fPs*1e3:.5f} mm = {h*fPs*1e6:.3f} ueV")
gapH = relH * fH  # C1 without density, measured g: change of the F = 1 - F = 0 splitting
out.append(f"H  C1 gap: {gapH/1e6:+.4f} MHz ({h*gapH*1e6:+.6f} ueV) of the splitting;"
           f" F = 1 {h*gapH/4*1e6:+.6f} ueV, F = 0 {-3*h*gapH/4*1e6:+.6f} ueV")
mag = float(((1 + A_E) ** 2)) * alpha**4 * mc2 / 3 / h   # magnetic T - S of Ps 1S, Hz
gapP = relP * mag
out.append(f"Ps magnetic T - S (g measured) {mag/1e9:.4f} GHz; C1 gap {gapP/1e6:+.2f} MHz ({h*gapP*1e6:+.4f} ueV) of T - S:"
           f" o-Ps (triplet) {h*gapP/4*1e6:+.4f} ueV, p-Ps (singlet) {-3*h*gapP/4*1e6:+.4f} ueV"
           f"  -- {gapP/fPs*1e6:+.0f} ppm of the measured HFS")
text = "\n".join(out)
print(text)
open("out/hfs682.txt", "w").write(text + "\n")

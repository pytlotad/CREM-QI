# Audit 756: analysis of the p585 emission ledger -> out/ana756.txt
import os, re, glob, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
hb = 1.054571817e-34; eV = 1.602176634e-19; Ry = 13.605693122994; hc = 1239.84198
Bp = {'H': Ry * eV / (1 + 1 / 1836.15267343), 'Ps': Ry * eV * 0.5}
pat = re.compile(r'\| B E (\S+) L (\S+).*\| A E (\S+) L (\S+)')
okA = okB = okC = True; rows = []
for f in sorted(glob.glob(os.path.join(O, '*.err'))):
    tag = os.path.basename(f)[:-4]; pair, ns, Ls, W = tag.split('_'); B = Bp[pair]
    end = [l for l in open(f) if l.startswith('END')]; ph = []
    for l in open(f):
        if l.startswith('CREM_LEDGER t'):
            m = pat.search(l); Eb, Lb, Ea, La = map(float, m.groups()); p = l.split()
            ph.append((float(p[2]), float(p[4]) / eV, (B / -Eb) ** .5, Lb / hb, (B / -Ea) ** .5, La / hb))
    pr(f"{tag}: {end[0].strip() if end else 'brak END'}; fotonow w ksiedze {len(ph)}")
    for t, Eg, nb, lb, na, la in ph:
        fr = abs((nb - round(nb)) - (na - round(na))) if True else 0
        dfr = abs((nb % 1) - (na % 1)); dfr = min(dfr, 1 - dfr); dL = min(abs(la - (lb - 1)), abs(la - (lb + 1)))
        lock_n = abs(nb - round(nb)); l = round(lb - 0.5); lock_L = abs(lb - (l + 0.5))
        pr(f"   t {t:.4e} s  E_gamma {Eg:.6f} eV ({hc/Eg:.4f} nm, {Eg/4.135667696e-15:.5e} Hz)  przed n {nb:.6f} L {lb:.6f}  po n {na:.6f} L {la:.6f}"
           f"  |d frac n| {dfr:.1e} |dL-+1| {dL:.1e}")
        if W == 'W0': okA &= dfr < 1e-4 and dL < 1e-3
        else:
            okB &= lock_n < 1e-3 and lock_L < 1e-2
            if abs(na - 1) < 1e-3 and abs(nb - 2) < 1e-3:
                ref = 0.75 * B / eV; okC &= abs(Eg - ref) < 1e-4
                pr(f"     linia 2 -> 1: {Eg:.6f} eV wobec Bohra {ref:.6f} eV (roznica {Eg-ref:+.2e} eV)")
    rows.append((tag, len(ph)))
pr(f"H756a (W0: foton zachowuje czesci ulamkowe n i L): {'TAK' if okA else 'NIE'}")
pr(f"H756b (W1: przed fotonem n calkowite, L = l + 1/2): {'TAK' if okB else 'NIE'}")
pr(f"H756c (W1: linia 2 -> 1 = Bohr do 1e-4 eV): {'TAK' if okC else 'NIE'}")
pr("Pomiar: Lyman alfa wodoru 10,19881 eV (121,5668 nm; NIST, srodek ciezkosci ze struktura subtelna); Ps 2P -> 1S 5,1021 eV (243,00 nm)")

# Audit 768: replace the rho3 table of precessionSpinOrbitEnergy (modules/crem_collapse.hpp) by the k = 2 table of audit 767.
import os, re
D = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(D, '..', '..', '..')
T = {}
for x in open(os.path.join(D, '..', 'section-767', 'out', 'rho3_k2_fine2_table.txt')):
    if x[0] != '#':
        p = x.split(); T[int(p[0]), int(p[1])] = p[2]
assert len(T) == 28
row = lambda n: '{' + ', '.join(T.get((n, l), '0.0') for l in range(8)) + '}'
new = 'static const double rho3[9][8]={\n' + ',\n'.join('        ' + row(n) for n in range(9)) + '};'
H = os.path.join(R, 'modules', 'crem_collapse.hpp'); s = open(H, encoding='utf8').read()
s2, c = re.subn(r'static const double rho3\[9\]\[8\]=\{.*?\}\};', lambda m: new, s, count=1, flags=re.S); assert c == 1
a = """// rho3 from the EBK wave of principle I at order hbar^2 (audit 737, one Schwarzian iteration; table from
// audits/data/section-763/gen763.py, n = 2..8; rho3/QM - 1 between -5.3e-5 and -1e-7).  Outside the table the caller keeps
// Sommerfeld.  Offline (762): 2P3/2 - 2P1/2 = 10967.97 MHz vs measured 10969.045.
"""
b = """// rho3 from the EBK wave of principle I at order hbar^4 (audits 764-767, two Schwarzian iterations, window method; table from
// audits/data/section-767/out/rho3_k2_fine2_table.txt via section-768/gen768.py, n = 2..8; |rho3/QM - 1| <= 2.5e-6).  Outside the
// table the caller keeps Sommerfeld.  Offline (767): 2P3/2 - 2P1/2 = 10968.715 MHz vs measured 10969.045.
"""
assert s2.count(a) == 1; open(H, 'w', encoding='utf8').write(s2.replace(a, b))

# Audit 766: verdicts from the logs of run 3 (2p, 5p, 8p, 3d; OOM-killed in 8d) + dokoncz766.py (8d) -> out/ana766.txt
import os, re
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out')
L3 = [x.rstrip('\n') for x in open(os.path.join(O, 'run3_2p_5p_8p_3d.log'))]; L8 = [x.rstrip('\n') for x in open(os.path.join(O, 'run_8d.log'))]
lines = [x for x in L3 if not x.startswith('  8d')] + L8
g = lambda pat, s: float(re.search(pat, s).group(1))
V = {}
for s in ('2p', '5p', '8p', '3d', '8d'):
    f1 = next(x for x in lines if x.startswith(f'  {s} fine 1')); f2 = next(x for x in lines if x.startswith(f'  {s} fine 2'))
    rm = next(x for x in lines if x.startswith(f'  {s} rmin')); gl = next(x for x in lines if x.startswith(f'  {s} gl'))
    k = r'\|K3-k1\| (\S+)'
    V[s] = dict(K1=g(k, f1), K2=g(k, f2), k21=g(r'k2 (\S+)', f1), k22=g(r'k2 (\S+)', f2), a6=g(r'rmin 1e-6: k2 (\S+)', rm), a4=g(r'rmin 1e-4: k2 (\S+)', rm), g2=g(r'k2 (\S+) \(gl', gl))
a = all(v['K2'] < 1e-8 for v in V.values()) and all(V[s]['K1'] >= 3 * V[s]['K2'] for s in ('5p', '8p'))
ps = ('2p', '5p', '8p'); ds = ('3d', '8d')
dm = {s: max(abs(V[s]['k22'] - V[s]['k21']), abs(V[s]['a6'] - V[s]['k21']), abs(V[s]['a4'] - V[s]['k21'])) for s in ps}
b = all(x < 3e-7 for x in dm.values())
gp = {s: abs(V[s]['g2'] - V[s]['k21']) for s in ps}; gd = {s: abs(V[s]['g2'] - V[s]['k21']) for s in ds}
c = all(x > 1e-6 for x in gp.values()); cd = all(x < 3e-7 for x in gd.values())
lines += [f"H766a (K3 przy fine 2 < 1e-8, spadek >= 3x dla 5p, 8p): {'TAK' if a else 'NIE'}; " + ' '.join(f"{s} {V[s]['K1']:.1e}->{V[s]['K2']:.1e}" for s in V),
          f"H766b (resztka p niezalezna od siatki i odciecia do 3e-7): {'TAK' if b else 'NIE'}; najwieksze zmiany " + ' '.join(f"{s} {x:.1e}" for s, x in dm.items()),
          f"H766c (|gl - in| > 1e-6 dla p; d < 3e-7): p {'TAK' if c else 'NIE'} (" + ' '.join(f"{s} {x:.1e}" for s, x in gp.items()) + f"), d {'TAK' if cd else 'NIE'} (" + ' '.join(f"{s} {x:.1e}" for s, x in gd.items()) + ")"]
open(os.path.join(O, 'ana766.txt'), 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines[-3:]))

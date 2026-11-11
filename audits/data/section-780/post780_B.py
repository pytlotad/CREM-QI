# Audit 780, POST HOC (nie zmienia werdyktow): odleglosc od EBK wg zwyciezcy -- 774G w 3 tau i na koncu (6 tau) wobec B w 3 tau; C wobec 775B05 -> out/post780_B.txt
import os, sys, numpy as np, io, contextlib
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
with contextlib.redirect_stdout(io.StringIO()): import ana780 as a
R = a.R; ids = a.ids
q4 = {i: a.state(a.load(os.path.join(D, '..', 'section-774', 'out', 'runs', f"G_{i}.txt")), 20000) for i in ids}
print("id  zwyc.(B)  d 774G(3tau)  d 774G(6tau)  d B(3tau)   | zwyc.(775B05) d B05(3tau) ; zwyc.(C) d C(3tau)")
for i in ids:
    b = R['B (GT 0,0001)'][i]; g = R['774G (GT 0,0003)'][i]; c = R['C (beta 0,0125, GT 0,000075)'][i]; p = R['775B05 (beta 0,05, GT 0,0003)'][i]
    print(f"{i}  {str(b['dom']):9s} {g['d']:.4f}{'*' if g['one'] else ' '}       {q4[i]['d']:.4f}{'*' if q4[i]['one'] else ' '}       {b['d']:.4f}{'*' if b['one'] else ' '}    | {str(p['dom']):9s} {p['d']:.4f} ; {str(c['dom']):9s} {c['d']:.4f} w {c['w']:.2f}")
O = [i for i in ids if R['B (GT 0,0001)'][i]['one'] and q4[i]['one']]
print(f"pary jednomodowe w obu (774G 6 tau, B 3 tau): {len(O)}; mediana d 774G {np.median([q4[i]['d'] for i in O]):.4f}, B {np.median([R['B (GT 0,0001)'][i]['d'] for i in O]):.4f}; B mniejsze w {sum(R['B (GT 0,0001)'][i]['d'] < q4[i]['d'] for i in O)}/{len(O)}")
nr0 = lambda arm: sum(R[arm][i]['dom'][0] == 0 for i in ids)
print(f"udzial n_r = 0 wsrod zwyciezcow: 775B05 {nr0('775B05 (beta 0,05, GT 0,0003)')}/20, C {nr0('C (beta 0,0125, GT 0,000075)')}/20; 774G {nr0('774G (GT 0,0003)')}/20, B {nr0('B (GT 0,0001)')}/20")

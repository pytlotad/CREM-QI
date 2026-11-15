# Audit 795 POST HOC (written after ana795 results): (1) which t795 lines miss H795t and the size of the circular-torus residual against
# the two-body scale nu alpha^2/n^2 (nu = m1 m2/M^2); (2) what the production start changed in the regression sets (old vs new out.txt).
# usage: python3 post795.py   (R795REG = scratchpad r795reg with reg/{old,new}/<set>/out.txt)
import os
D = os.path.dirname(os.path.abspath(__file__))
REG = os.environ.get('R795REG', '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r795reg') + '/reg'
A = 7.2973525693e-3; ME = 9.1093837015e-31; MP = 1.67262192369e-27
NU = {'positron,electron': 0.25, 'proton,electron': MP*ME/(MP + ME)**2}
pr = lambda *a: print(*a, flush=True)
pr('=== t795: linie poza progiem H795t (|E_orb/E_S - 1| <= 1e-9) ===')
worst = {}
for ln in open(os.path.join(D, 'out', 't795.txt')):
    w = ln.split(); e = float(w[w.index('nowy:') + 2]); n = int(w[w.index('n') + 1]); L = float(w[w.index('L') + 1]); dip = w[2]
    key = (w[0], n, L); worst[key] = max(worst.get(key, 0.0), abs(e))
    if abs(e) > 1e-9:
        x = NU[w[0]]*A*A/n**2
        pr(f'  {w[0]} dip {dip} n {n} L {L:g} (kolowy: {L == n}): E_orb/E_S - 1 {e:+.4e} = {e/x:+.4f} x nu alfa^2/n^2 ({x:.4e})')
pr('  maks. |E_orb/E_S - 1| poza kolowymi: %.1e' % max(v for k, v in worst.items() if k[2] != k[1]))
pr('=== regresja: old vs new (out.txt), linie rozne ===')
for s in ('1', '2', 'L2', '5Ps', '5H', '6'):
    a = open(f'{REG}/old/{s}/out.txt').read().split('\n'); b = open(f'{REG}/new/{s}/out.txt').read().split('\n')
    d = [(x, y) for x, y in zip(a, b) if x != y]
    pr(f'  zestaw {s}: {len(d)} linii rozne (dlugosc {len(a)} / {len(b)})')
    for x, y in d: pr(f'    old: {x.strip()[:150]}\n    new: {y.strip()[:150]}')
    for ln in b:
        if ln.startswith('Annihilation') or 'E[T] from the preparation' in ln: pr(f'    (bez zmian) {ln.strip()[:170]}')

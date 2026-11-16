# Audit 797 POST HOC (written after ana797 results): do the run-to-run differences of Omega (seeds, start 791 -> 797, F = 0 vs F = 1, 2048 vs 512)
# take discrete per-torus values (a whole periapsis kick in or out of the fixed-time mechB0 window) rather than a continuous spread?
# usage: python3 post797.py   (R791 as in ana797.py)
import sys, os, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
sys.path.insert(0, os.path.join(D, '..', 'section-699')); import ana699
R791 = os.environ.get('R791', '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r791x/runs')
TORI = [(3, 1), (3, 1.5), (3, 2), (3, 2.5), (3, 3), (2, 1), (2, 1.5), (2, 2)]
lab = lambda L: '%g' % L
A, C, A0, ME, MP = 7.2973525693e-3, 299792458.0, 5.29177210903e-11, 9.1093837015e-31, 1.67262192369e-27
MU = ME*MP/(ME + MP)
pr = lambda *a: print(*a, flush=True)
def om(path):
    m = ana699.parse(path)['inc']['mechB0']
    return np.linalg.norm(m['s21'] - m['s20'])/(m['dt']*np.linalg.norm(np.cross(m['L0'], m['s20'])))/(2*math.pi), m['dt']
for k in TORI:
    n, L = k; T = 2*math.pi*n**3*A0*ME/MU/(A*C)
    V = {}
    for s in '456':
        for t in ('proton', 'antiproton'):
            V[f'797 {t[0]}{s}'] = om(f'{O}/G/{t}_n{n}L{lab(L)}_{s}.err'); V[f'791 {t[0]}{s}'] = om(f'{R791}/G/{t}_n{n}L{lab(L)}_{s}.err')
    for t in ('proton', 'antiproton'): V[f'F0 {t[0]}5'] = om(f'{O}/G0/{t}_n{n}L{lab(L)}_5.err')
    if k == (3, 1):
        for t in ('proton', 'antiproton'): V[f'2048 {t[0]}5'] = om(f'{O}/Z/{t}_n3L1_5.err')
    lo = min(v[0] for v in V.values()); dt = list(V.values())[0][1]; N = dt/T
    pr(f'n{n}L{lab(L)}: okno dt / T_n = {N:.3f} obiegu, 1/N = {1/N:.2e}; Omega/min - 1 posortowane:')
    pr('   ' + ', '.join(f'{name} {v[0]/lo - 1:.2e}' for name, v in sorted(V.items(), key=lambda x: x[1][0])))
pr('=== zmiany pojedynczego biegu w jednostkach 1/N (N = 200 obiegow w oknie): 797/791 - 1 (start), F0/F1 - 1 (znak momentu), 2048/512 - 1 ===')
for k in TORI:
    n, L = k; cells = []
    for s in '456':
        for t in ('proton', 'antiproton'):
            cells.append(f'{t[0]}{s} {(om(f"{O}/G/{t}_n{n}L{lab(L)}_{s}.err")[0]/om(f"{R791}/G/{t}_n{n}L{lab(L)}_{s}.err")[0] - 1)*200:+.2f}')
    for t in ('proton', 'antiproton'):
        cells.append(f'F0/F1 {t[0]}5 {(om(f"{O}/G0/{t}_n{n}L{lab(L)}_5.err")[0]/om(f"{O}/G/{t}_n{n}L{lab(L)}_5.err")[0] - 1)*200:+.2f}')
    if k == (3, 1):
        for t in ('proton', 'antiproton'): cells.append(f'2048 {t[0]}5 {(om(f"{O}/Z/{t}_n3L1_5.err")[0]/om(f"{O}/G/{t}_n3L1_5.err")[0] - 1)*200:+.2f}')
    pr(f'  n{n}L{lab(L)}: ' + ', '.join(cells))

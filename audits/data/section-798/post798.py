# Audit 798 POST HOC (written after ana798 results): F = 0 / F = 1 and old / new start in units of 1/N (N = 200) per torus, with the apoapsis window;
# and the orbit-count estimate: a relative energy offset dE/E shifts the Kepler period by 1.5 dE/E, i.e. by 300 dE/E orbits over 200 orbits.
import sys, os, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
sys.path.insert(0, os.path.join(D, '..', 'section-699')); import ana699
TORI = [(3, 1), (3, 1.5), (3, 2), (3, 2.5), (3, 3), (2, 1), (2, 1.5), (2, 2)]
lab = lambda L: '%g' % L
pr = lambda *a: print(*a, flush=True)
def om(path):
    m = ana699.parse(path)['inc']['mechB0']
    return np.linalg.norm(m['s21'] - m['s20'])/(m['dt']*np.linalg.norm(np.cross(m['L0'], m['s20'])))/(2*math.pi)
pr('F = 0 / F = 1 - 1, x 200 (okno w apocentrum, ziarno 5):')
for t in ('proton', 'antiproton'):
    pr(f'  {t}: ' + ', '.join(f'n{k[0]}L{lab(k[1])} {(om(f"{O}/A0/{t}_n{k[0]}L{lab(k[1])}_5.err")/om(f"{O}/A/{t}_n{k[0]}L{lab(k[1])}_5.err") - 1)*200:+.3f}' for k in TORI))
pr('luka proton / antyproton - 1, x 200 (okno w apocentrum, F = 1):')
for s in '456':
    pr(f'  ziarno {s}: ' + ', '.join(f'n{k[0]}L{lab(k[1])} {(om(f"{O}/A/proton_n{k[0]}L{lab(k[1])}_{s}.err")/om(f"{O}/A/antiproton_n{k[0]}L{lab(k[1])}_{s}.err") - 1)*200:+.2f}' for k in TORI))
pr('szacunek liczby obiegow: przesuniecie okresu 1,5 dE/E, po 200 obiegach 300 dE/E obiegu --')
pr('  stary start n3L1 dE/E = 0,51 %% energii wiazania (794) -> %.2f obiegu; dryf energii protonu n3L1 -4,8e-3 w 200 obiegach (793, liniowo) -> %.2f obiegu' % (300*0.0051, 300*0.0048/2))

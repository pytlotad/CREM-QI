# Audit 793 POST HOC: seed-mean rule, Omega n^3 L^2 vs eccentricity, full-force channel gap at 512 / 2048.  usage: python3 post793.py
import io, contextlib, math
import numpy as np
with contextlib.redirect_stdout(io.StringIO()): import ana793 as A
pr = A.pr; O = A.O
W = {}
for n, L in A.TORI:
    W[(n, L)] = np.mean([A.omega(f'{O}/G/{t}_n{n}L{A.lab(L)}_{s}.err')[0] for t in ('proton', 'antiproton') for s in '456'])
pr('=== (a) regula ze srednich 6 biegow (CREM_NO_DIPOLE_FORCE) ===')
for n, l, st in ((2, 1, '2p'), (3, 1, '3p'), (3, 2, '3d')):
    od = l*W[(n, l)] - (l + 1)*W[(n, l + 1)]; q = od/W[(n, l + 0.5)]; qt = (l + 0.5)**2/(l*(l + 1)); rf = A.FS*16/(n**3*l*(l + 1))
    pr(f'  {st}: {A.units(od)} | {q/qt:.5f} x regula | {od/rf:.5f} x skalowanie pomiaru 2P ({(od - rf)/1e6:+.2f} MHz)')
pr('=== (b) Omega n^3 L^2 wobec mimosrodu e = sqrt(1 - L^2/n^2) (srednie 6 biegow; wzgledem n = 3 L = 3) ===')
ref = W[(3, 3)]*27*9
for n, L in sorted(A.TORI, key=lambda k: -math.sqrt(1 - k[1]**2/k[0]**2)):
    pr(f'  n {n} L {L:3}: e {math.sqrt(1 - L**2/n**2):.3f} | rp/a {1 - math.sqrt(1 - L**2/n**2):.3f} | Omega n^3 L^2 {W[(n, L)]*n**3*L**2/ref - 1:+.4%}')
pr('=== (c) luka kanalow n3L1 ziarno 5 z silami od momentow: 512 / 2048 ===')
for w in ('E512', 'E2048', 'E0'):
    a = A.omega(f'{O}/E/proton_{w}.err')[0]; b = A.omega(f'{O}/E/antiproton_{w}.err')[0]
    pr(f'  {w}: proton {a/1e6:.2f} MHz, antyproton {b/1e6:.2f} MHz -> {a/b - 1:+.4f}')

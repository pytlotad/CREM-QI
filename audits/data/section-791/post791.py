# Audit 791 POST HOC (written after the scored ana791 output): (a) per-torus mean and spread of Omega over the 6 runs (proton/antiproton x seeds 4-6)
# and the difference rule from the means; (b) the L = 1 window / step diagnostic (out/runs/D, runs/diag.sh), if present.
import os, math, glob
import numpy as np
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()): from ana791 import omega, units, TORI, lab, O, FS, pr   # silence ana791's scoring printout
W = {}
for n, L in TORI:
    v = [omega(f'{O}/G/{t}_n{n}L{lab(L)}_{s}.err') for t in ('proton', 'antiproton') for s in '456']
    om = np.array([a for a, _ in v]); ln = v[0][1]; W[(n, L)] = (om.mean(), ln)
    pr(f'  n {n} L {L:3}: srednia {units(om.mean())} | rozrzut (max/min-1) {(om.max()/om.min()-1)*100:.2f} % | sd/srednia {om.std(ddof=1)/om.mean()*100:.2f} % | Omega n^3 L^2 {om.mean()*n**3*ln**2/1e6:.3f} MHz')
k = [om*n**3*ln**2 for (n, L), (om, ln) in W.items()]; pr(f'  rozrzut Omega n^3 L^2 ze srednich: {(max(k)/min(k)-1)*100:.3f} %')
for n, l, st in ((2, 1, '2p'), (3, 1, '3p'), (3, 2, '3d')):
    F = lambda L: W[(n, L)][1]*W[(n, L)][0]
    od = F(l) - F(l + 1); oc = W[(n, l + 0.5)][0]; q = od/oc; qt = (l + 0.5)**2/(l*(l + 1)); rf = FS*16/(n**3*l*(l + 1))
    pr(f'  {st} ze srednich: {units(od)} | stosunek {q:.5f} wobec {qt:.5f} -> {q/qt:.5f} | wobec skalowania pomiaru {rf/1e6:.2f} MHz -> {od/rf:.5f} ({(od-rf)/1e6:+.2f} MHz)')
pr('=== diagnostyka L = 1 (n = 3, ziarno 5) ===')
for t in ('proton', 'antiproton'):
    base = omega(f'{O}/G/{t}_n3L1_5.err')[0]; pr(f'  {t} okno 200, krok 512: {base/1e6:.2f} MHz')
    for w in ('T400', 'S2048'):
        p = f'{O}/D/{t}_n3L1_5_{w}.err'
        if os.path.exists(p) and 'mechB0' in open(p).read(): a = omega(p)[0]; pr(f'  {t} {w}: {a/1e6:.2f} MHz ({a/base-1:+.4f} wobec bazy)')
# (c) J-frame check (POST HOC; refuted): the orbit recoils by the spin change (L1 - L0 = -(s21 - s20)/2, J = L + S_e conserved), so the spin
# precesses about J. But J x s = Ln0 L0 x s + s x s/2 = Ln0 L0 x s, so the 699 estimator already measures xi L; the recoil does not bias it.
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-699')); import ana699
pr('=== (c) odrzut orbity i estymator w ukladzie J (POST HOC, hipoteza obalona) ===')
dmax = 0.0
for f in sorted(glob.glob(f'{O}/G/*.err')):
    m = ana699.parse(f)['inc']['mechB0']; J = m['Ln0']*m['L0'] + 0.5*m['s20']
    a = np.linalg.norm(np.cross(J, m['s20'])); b = m['Ln0']*np.linalg.norm(np.cross(m['L0'], m['s20'])); dmax = max(dmax, abs(a/b - 1))
m = ana699.parse(f'{O}/G/proton_n3L1_5.err')['inc']['mechB0']
pr(f"  proton n3L1 ziarno 5: Ln0 (L1 - L0) = {np.round(m['Ln0']*(m['L1'] - m['L0']), 5)}, -(s21 - s20)/2 = {np.round(-(m['s21'] - m['s20'])/2, 5)}")
pr(f'  |J x s| / (Ln0 |L0 x s|) - 1 po 48 biegach G: max {dmax:.1e} -> estymator J = estymator 699')

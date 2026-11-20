# Audit 817 (= ana816 + H817c, cdarw817): is delta_6 of 815 a property of the dynamics or the start label?  Action-angle variables of the Darwin Lagrangian (cdarw817)
# on the 815 starts; written before computing.  Usage: python3 ana817.py -> out/ana817.txt
import os, sys, math, tarfile
import numpy as np, mpmath as mp
D = os.path.dirname(os.path.abspath(__file__)); sys.path[:0] = [D, os.path.join(D, '..', 'section-815'), os.path.join(D, '..', 'section-814')]
from cdarw817 import Pair, HBAR
import pred815 as P15
from pred814 import parse_start
pr = lambda *a: print(*a, flush=True)
H = 2*mp.pi*HBAR; HEV, C = 4.135667696e-15, 299792458.0; al = 1/137.035999177; NU4 = 5474.642e6
units = lambda f: f'{f/1e6:.4f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
x, w = np.polynomial.legendre.leggauss(5); u = 0.75 + 0.25*x[::-1]
NODES = dict(zip([f'u5{c}' for c in 'abcde'], zip(1/u, 0.25*w[::-1])))
nu_of = lambda Om: sum(w*L**2*Om[k] for k, (L, w) in NODES.items())/(2*math.pi)
T815 = os.path.join(D, '..', 'section-815', 'runs', 'runs815.tar.xz')
def starts(prefix):
    out = {}
    with tarfile.open(T815) as t:
        for k in NODES:
            m = [q for q in t.getnames() if q.endswith(prefix.format(k))][0]
            tmp = os.path.join(D, 'out', f'_{k}.err'); open(tmp, 'wb').write(t.extractfile(m).read())
            out[k] = parse_start(tmp); os.remove(tmp)
    return out
SP = {'o-Ps (PHEN=2)': starts('S025/F1_ops_{}_s8192_5.err'), 'p-Ps (PHEN=1)': starts('PP025/F0_pps_{}_s2048_5.err')}
s0 = SP['o-Ps (PHEN=2)']['u5a']; m1, m2, kq, c = s0['m1'], s0['m2'], s0['kq'], s0['c']
pr(f'stale (PAIR CONS814): m1 {m1:.10e}, m2 {m2:.10e}, kq1q2 {kq:.10e}, c {c:.1f}')
# ---------------- H817t: no Darwin, exact closed form
P0 = Pair(m1, m2, kq, c, dar=0); k = -P0.kq; dt = 0
for kk, st in SP['o-Ps (PHEN=2)'].items():
    E, Lc = P0.start_invariants(st)[:2]
    A_ = E**2/(4*P0.c**2) - P0.m1**2*P0.c**2; B = E*k/(4*P0.c**2)
    Jc = B/mp.sqrt(-A_) - mp.sqrt(Lc**2 - k**2/(4*P0.c**2)); Jn = P0.Jr(E, Lc); dt = max(dt, abs(Jn - Jc)/HBAR)
    pr(f'  {kk}: bez Darwina J_r/hbar kwadratura {mp.nstr(Jn/HBAR, 20)}, wzor {mp.nstr(Jc/HBAR, 20)}')
pr(f"H817t (A = 0: kwadratura = wzor zamkniety do 1e-20 hbar): {'TAK' if dt <= 1e-20 else 'NIE'} (max {mp.nstr(dt, 3)})")
Ekin = lambda Jr, L: 2*P0.m1*P0.c**2/mp.sqrt(1 + k**2/(4*P0.c**2*(Jr + mp.sqrt(L**2 - k**2/(4*P0.c**2)))**2))
nu_kin = float((Ekin(0, 2*HBAR) - Ekin(HBAR, HBAR))/H)
ec = P0.ecirc(2*HBAR); ecf = Ekin(0, 2*HBAR); dcirc = abs((ec - P0.M*P0.c**2)/(ecf - P0.M*P0.c**2) - 1)
pr(f"H817c (A = 0: E_circ(2 hbar) = wzor zamkniety do 1e-12): {'TAK' if dcirc <= 1e-12 else 'NIE'} ({mp.nstr(dcirc, 3)})")
# ---------------- H817s: one-body limit
PS = Pair(m1, 1e12*m1, kq, c, dar=1)
nu_s = (PS.etorus(2, 2*HBAR) - PS.etorus(2, HBAR))/H
Es = lambda Jr, L: P0.m1*P0.c**2/mp.sqrt(1 + (k/(P0.c*(Jr + mp.sqrt(L**2 - k**2/P0.c**2))))**2)
nu_sc = (Es(0, 2*HBAR) - Es(HBAR, HBAR))/H; ds = abs(nu_s/nu_sc - 1)
pr(f"H817s (m2 = 1e12 m1: roznica torusow = Sommerfeld do 1e-9): {'TAK' if ds <= 1e-9 else 'NIE'} ({mp.nstr(nu_s/1e6, 15)} vs "
   f"{mp.nstr(nu_sc/1e6, 15)} MHz, {mp.nstr(ds, 3)})")
# ---------------- per species
okQ = okD = okL = True; RES = {}
for nm, SS in SP.items():
    st = SS['u5a']; P = Pair(st['m1'], st['m2'], st['kq'], st['c'], dar=1)
    E1 = P.etorus(2, HBAR); E2 = P.etorus(2, 2*HBAR); nu_true = float((E2 - E1)/H)
    Omq = {}
    for kk, (L, ww) in NODES.items():
        El = P.etorus(2, L*HBAR, E0=None); Omq[kk] = float(P.omegas(El, L*HBAR)[0])
    nu_q = nu_of(Omq); okQ &= abs(nu_q - nu_true) <= 1e-4*1e6
    pr(f'{nm}: nu_true = {units(nu_true)}; kwadratura 808 na torusach prawdziwych {nu_q/1e6:.6f} MHz ({(nu_q - nu_true)/1e6:+.2e})')
    Oe, Od = {}, {}; neff = 0
    for kk, st in SS.items():
        E, Lc = P.start_invariants(st)[:2]; Od[kk] = float(P.omegas(E, Lc)[0]); Oe[kk] = P15.omega(st, 0.0)[0]
        ne = (P.Jr(E, Lc) + Lc)/HBAR - 2; neff = max(neff, abs(ne)); Llab = NODES[kk][0]
        rd = Oe[kk]/Od[kk] - 1; okD &= abs(rd) <= 5e-6
        pr(f'  {kk}: Omega_silnik(eps -> 0) {Oe[kk]:.10e}, Omega_Darwin(start) {Od[kk]:.10e} ({rd:+.2e}); (J_r + L_can)/hbar - 2 = '
           f'{mp.nstr(ne, 4)}; L_can/L_etykieta - 1 = {mp.nstr(Lc/(Llab*HBAR) - 1, 4)}')
    nu815 = nu_of(Oe); shift = nu815 - nu_true; okL &= neff >= 1e-6 and abs(shift) >= 0.010e6
    RES[nm] = (nu_true, nu815, shift)
    pr(f'  nu_815(eps -> 0) = {nu815/1e6:.4f} MHz; przesuniecie etykiety nu_815 - nu_true = {shift/1e6:+.4f} MHz; max |n_eff - 2| = {mp.nstr(neff, 3)}')
    pr(f'  RAPORT: delta_6,true = {(nu_true - NU4)/1e6:+.4f} MHz (815: {(nu815 - NU4)/1e6:+.4f}; jednocialowy +0,182); rho = '
       f'{nu_true/NU4:.7f} = 1 {(nu_true/NU4 - 1)/al**2:+.3f} alfa^2; kinetyka bez Darwina {nu_kin/1e6:.4f} MHz, czesc Darwina '
       f'{(nu_true - nu_kin)/1e6:.4f} MHz')
pr(f"H817q (metoda 808 na torusach prawdziwych = nu_true do 1e-4 MHz): {'TAK' if okQ else 'NIE'}")
pr(f"H817d (GLOWNA 1: Omega_silnik = Omega_Darwin(start) do 5e-6, o-Ps i p-Ps): {'TAK' if okD else 'NIE'}")
pr(f"H817l (GLOWNA 2: starty poza torusami n = 2 >= 1e-6 i przesuniecie etykiety >= 0,010 MHz): {'TAK' if okL else 'NIE'}")
val = dt <= 1e-20 and dcirc <= 1e-12 and ds <= 1e-9 and okQ
pr('DECYZJA: ' + ('wynik niewazny (walidacja NIE)' if not val else
   ('delta_6 orbity produkcyjnej = lagranzjan Darwina; wartosc modelu delta_6,true (torusy prawdziwe), 815 zawiera etykiete startu' if okD and okL else
    'delta_6 orbity produkcyjnej = lagranzjan Darwina; etykieta startu bez znaczenia ponad progi' if okD else
    'silnik rozni sie od lagranzjanu Darwina w rzedzie alfa^6 tempa apsyd (przyspieszenia zredukowane); raport roznicy')))

# Audit 811: verdicts H811f, H811r, H811c (part A) and H811d, H811i, H811e (part B) -> out/ana811.txt
# Inputs: out/p811_new.txt, out/p811_off.txt (runs/all.sh), ../section-810/out/p810_{new,new_off}.txt
import os, math
import sympy as sp
from mpmath import mp, mpf, sqrt as msqrt, quad, cos, pi as mpi
mp.dps = 40
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out')
h = 6.62607015e-34; HEV = 4.135667696e-15; C = 299792458.0
pr = lambda *a: print(*a, flush=True)
rd = lambda f: open(f).read().splitlines()
# ---------------- part A
al = mpf(1)/mpf('137.035999177'); me = mpf('9.1093837139e-31'); mpp = mpf('1.67262192595e-27'); cc = mpf(299792458); hh = mpf('6.62607015e-34')
mu = me*mpp/(me + mpp); lam = mpf('2.00231930436092') - 1
def ES(n, L): return mu*cc**2*(1/msqrt(1 + al**2/(n - L + msqrt(L**2 - al**2))**2) - 1)
def E4(n, L): return -mu*cc**2*al**2/(2*n**2) - mu*cc**2*al**4/(2*n**4)*(n/L - mpf(3)/4)
def R(n, L): return float((ES(n, L) - E4(n, L))/hh/10**6)
def lev(n, l, j2):
    if l == 0: return R(n, 1)
    k = (j2 + 1)//2; return R(n, l + 0.5) + float(lam)*(R(n, k) - R(n, l + 0.5))
def levJ(f):
    out = {}
    for w in (x.split() for x in rd(os.path.join(O, f)) if x.startswith('levelJ ')):
        out[(w[1], int(w[9]))] = (int(w[3]), int(w[5]), int(w[7]), float(w[13]), float(w[16]))
    return out
N, F_ = levJ('p811_new.txt'), levJ('p811_off.txt')
names = sorted({k[0] for k in N}, key=lambda s: (int(s[0]), 'SPDF'.index(s[1]), s[2:]))
dF = 0.0; dR = 0.0
pr('Czesc A: zmiana poziomu new - off [MHz] (z J, double)')
for nm in names:
    n, l, j2 = N[(nm, 1)][:3]
    d1 = (N[(nm, 1)][3] - F_[(nm, 1)][3])/h/1e6; d0 = (N[(nm, 0)][3] - F_[(nm, 0)][3])/h/1e6
    m = lev(n, l, j2); dF = max(dF, abs(d1 - d0)); dR = max(dR, abs(d1 - m), abs(d0 - m))
    pr(f'  {nm}: F = 1 {d1:+.12f}, F = 0 {d0:+.12f}, roznica {d1 - d0:+.2e}, model {m:+.12f}')
pr(f"H811f (|zmiana F=1 - zmiana F=0| <= 1e-9 MHz, {len(names)} poziomow): {'TAK' if dF <= 1e-9 else 'NIE'} (max {dF:.2e} MHz)")
pr(f"H811r (zmiana = model 810 do 1e-6 MHz, oba F): {'TAK' if dR <= 1e-6 else 'NIE'} (max {dR:.2e} MHz)")
def mhz(f, tag, emit):
    out = {}
    for w in (x.split() for x in rd(f) if x.startswith(tag)):
        out[(w[1], w[3])] = (w[emit], w[emit + 3])
    return out
okC = True
for a, b in (('p811_new.txt', 'p810_new.txt'), ('p811_off.txt', 'p810_new_off.txt')):
    A = mhz(os.path.join(O, a), 'levelMHz ', 5); B = mhz(os.path.join(D, '..', 'section-810', 'out', b), 'level ', 7)
    okC &= all(A[k] == B[k] for k in B) and len(B) == 18
pr(f"H811c (9 poziomow 810 w MHz %.9f znak w znak jak p810): {'TAK' if okC else 'NIE'}")
pr(f"DECYZJA A: {'zmiana 810 POTWIERDZONA jako domyslna' if dF <= 1e-9 and dR <= 1e-6 else 'zmiana 810 do opcji (opt-in) w 812'}")
# ---------------- part B
pr('Czesc B: <r^-3> Borna a roznica torusow orbity 1PN')
m_, k_, a_, L_, nu, E0, r = sp.symbols('mu k a L nu E0 r', positive=True)
A1, A2 = 1/a_, sp.sqrt(m_*k_)/(a_**sp.Rational(3, 2)*L_); A3 = m_*k_*A2/L_**2      # Kepler orbit averages <1/r>, <1/r^2>, <1/r^3>
E0v = -k_/(2*a_)
p4 = 4*m_**2*(E0v**2 + 2*E0v*k_*A1 + k_**2*A2)                                  # <p^4>, p^2 = 2 mu (E0 + k/r)
p2r = 2*m_*(E0v*A1 + k_*A2)                                                      # <p^2/r>
dar = 2*p2r - L_**2*A3                                                           # <(p^2 + p_r^2)/r>, p_r^2 = p^2 - L^2/r^2
c_ = sp.symbols('c', positive=True)
H1 = -(1 - 3*nu)*p4/(8*m_**3*c_**2) - k_*nu*dar/(2*m_**2*c_**2)
Ldep = sp.simplify(sp.expand(H1 - H1.subs(L_, 1)))                               # L-dependent part
target = sp.simplify(-(k_**2/(2*m_*c_**2))*(A2 - A2.subs(L_, 1)))
okDs = sp.simplify(Ldep - target) == 0 and sp.diff(sp.simplify(sp.diff(H1, L_)), nu) == 0
pr(f'  sympy: dH1/dL = {sp.simplify(sp.diff(H1, L_))};  d/dnu tego = {sp.simplify(sp.diff(sp.diff(H1, L_), nu))}')
# numeric check of the Kepler averages (eccentric anomaly quadrature)
okDn = True
for e in (mpf('0.3'), mpf('0.75')):
    a1 = mpf(1); Lq = msqrt(1 - e**2)                                            # mu = k = 1
    avg = lambda f: quad(lambda E: f(a1*(1 - e*cos(E)))*(1 - e*cos(E)), [0, 2*mpi])/(2*mpi)
    for f, ex in ((lambda rr: 1/rr, 1/a1), (lambda rr: rr**-2, 1/(a1**1.5*Lq)), (lambda rr: rr**-3, 1/(a1**1.5*Lq**3))):
        okDn &= abs(avg(f) - ex) < mpf('1e-30')
pr(f"H811d (czlon 1PN zalezny od L = -(k^2/(2 mu c^2)) <r^-2>, niezalezny od nu; srednie Keplera kwadratura do 1e-30): {'TAK' if okDs and okDn else 'NIE'}")
def born(mu_, k, hb, a, L): return mu_*k/(L*hb)*(msqrt(mu_*k)/(a**mpf(1.5)*(L - hb/2)) - msqrt(mu_*k)/(a**mpf(1.5)*(L + hb/2)))
def Eorb(mu_, k, cv, a, L): return -(k**2/(2*mu_*cv**2))*msqrt(mu_*k)/(a**mpf(1.5)*L)
hb = hh/(2*mpi); kk = al*hb*cc; worst = mpf(0)
for nm, m1, m2 in (('wodor', mpp, me), ('muonium', mpf('1.883531627e-28'), me), ('Ps', me, me)):
    mu_ = m1*m2/(m1 + m2); aB = hb**2/(mu_*kk)
    for l in range(1, 7):
        for n in range(l + 1, 9):
            a = n*n*aB; Ll = (l + mpf(1)/2)*hb
            rhs = 2*mu_**2*cc**2*(Eorb(mu_, kk, cc, a, (l + 1)*hb) - Eorb(mu_, kk, cc, a, l*hb))/(kk*hb*Ll/hb)/hb
            qm = 1/(aB**3*n**3*l*(l + mpf(1)/2)*(l + 1))
            worst = max(worst, abs(rhs/born(mu_, kk, hb, a, Ll) - 1), abs(rhs/qm - 1))
pr(f"H811i (Born = 2 mu^2 c^2 dE_orb/(k hbar L_l) = QM, l 1..6, n do 8, wodor / muonium / Ps, do 1e-30): {'TAK' if worst < mpf('1e-30') else 'NIE'} (max {float(worst):.1e})")
pairs = {}; we = 0.0
for w in (x.split() for x in rd(os.path.join(O, 'p811_new.txt'))):
    if w[0] == 'pair': pairs[w[1]] = [mpf(w[i]) for i in (3, 5, 7, 9, 11, 13)]
    if w[0] == 'born':
        m1, m2, k, hbe, ce, aP = pairs[w[1]]; mu_ = m1*m2/(m1 + m2); n, l = int(w[3]), int(w[5]); a = n*n*aP; Ll = (l + mpf(1)/2)*hbe
        rhs = 2*mu_**2*ce**2*(Eorb(mu_, k, ce, a, (l + 1)*hbe) - Eorb(mu_, k, ce, a, l*hbe))/(k*Ll*hbe)
        we = max(we, abs(float(mpf(w[7])/rhs - 1)))
pr(f"H811e (bornInverseCubeRadius silnika = roznica torusow, wodor i Ps, n 2..4): {'TAK' if we <= 1e-12 else 'NIE'} (max {we:.1e})")
m1, m2, k, hbe, ce, aP = pairs['electron,positron']; mu_ = m1*m2/(m1 + m2)
dE = Eorb(mu_, k, ce, 4*aP, 2*hbe) - Eorb(mu_, k, ce, 4*aP, hbe); f = float(dE/hh/1e6)
pr(f'PRZEWIDYWANIE 812 (Ps, orbita bez sil od momentow, J_r + L = 2): E_orb(2, 2) - E_orb(2, 1) = {f:.3f} MHz = {f*1e6*HEV*1e6:.4f} ueV = {C/(f*1e6)*100:.4f} cm '
   f'(= m_e c^2 alfa^4/64: {float(me*cc**2*al**4/64/hh/1e6):.3f} MHz); <r^-3>(2p Ps) = 2 mu^2 c^2 dE/(k hbar L_1) = 1/(24 a_Ps^3)')
m1, m2, k, hbe, ce, aP = pairs['proton,electron']; mu_ = m1*m2/(m1 + m2)
fH = float((Eorb(mu_, k, ce, 4*aP, 2*hbe) - Eorb(mu_, k, ce, 4*aP, hbe))/hh/1e6)
pr(f'  dla porownania wodor (rzad alfa^4): {fH:.3f} MHz; mechanika 808/809 z alfa^6: 10 943,7009 MHz')

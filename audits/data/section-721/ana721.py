"""Audit 721: hydrogen Lamb shift as an import (#25), from the CODATA 2006 theory of hydrogen levels (Mohr, Taylor, Newell, RMP 80, 633 (2008),
arXiv 0801.0028, Sec. IV.A.1, Eqs. (11), (12), (19)-(21), (22)-(25), (29)-(31), (51), (56); reduced-mass prescriptions as stated there).
Import per level (n, l), Z = 1, x = m_r/m_e, r = m_e/m_N, Lr = ln[(m_e/m_r) alpha^-2]:
  SE   (alpha/pi) alpha^4/n^3 x^3 [ (4/3) Lr d_l0 - (4/3) ln k0 + (10/9) d_l0 + (139/32 - 2 ln 2) pi alpha d_l0 ]
       (the magnetic-moment term -1/[2 kappa (2l+1)] of A40 is NOT imported: the engine has it since 516 via lambda = |g| - 1, import 21)
  VP   (alpha/pi) alpha^4/n^3 x^3 [ -(4/15) + (5 pi/48) alpha ] d_l0
  2L   (alpha/pi)^2 alpha^4/n^3 x^3 [ (3 pi^2/2) ln 2 - 10 pi^2/27 - 2179/648 - (9/4) zeta(3) ] d_l0   (l >= 1: magnetic line, in the engine's a_e)
  ES   x^3 r alpha^5/(pi n^3) [ (1/3) d_l0 ln alpha^-2 - (8/3) ln k0 - (1/9) d_l0 - (7/3) a_n - 2 d_l0 (mN^2 ln(me/mr) - me^2 ln(mN/mr))/(mN^2 - me^2) ]
  RR   x^3 r alpha^6/(pi^2 n^3) [ 6 zeta(3) - 2 pi^2 ln 2 + 35 pi^2/36 - 448/27 + (2/3) pi alpha ln^2 alpha^-2 ] d_l0
  NS   (2/3) x^3 alpha^2/n^3 (alpha R_p/lambdabar_C)^2 d_l0,  R_p = 0.84087(39) fm (Antognini et al., Science 339, 417 (2013), muonic hydrogen)
  all times m_e c^2.  Omitted: alpha (Z alpha)^6 and higher (A62, A61, G_SE, G_VP), B50 and higher two-loop, three-loop, D60 recoil, nuclear polarization.
Engine values (p721, before the change) from out/p721_before.txt.  usage: python3 ana721.py [p721 after-change output] -> pred721.txt / ana721.txt."""
import math, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
def rd(fn):
    d = {}
    for line in open(fn):
        k, v = line.split(); d[k] = float(v)
    return d
E = rd(os.path.join(HERE, 'out', 'p721_before.txt'))
a = 1/E['alpha_inv']; mc2 = E['me_c2_MHz']; mpme = E['mp_over_me']
x = mpme/(1 + mpme); r = 1/mpme
lnk0 = {}
for line in open(os.path.join(HERE, 'lnk0_n30.tsv')):
    n, l, v = line.split(); lnk0[(int(n), int(l))] = float(v)
def K(n, l): nn = min(n, 30); return lnk0[(nn, min(l, nn - 1))]
z3 = 1.2020569031595942
lamC = 386.15926796e-15   # m, reduced Compton wavelength of the electron (CODATA)
Rp = 0.84087e-15
def parts(n, l):
    s = 1.0 if l == 0 else 0.0; Lr = math.log(1/(x*a*a)); La = math.log(1/(a*a))
    pre = (a/math.pi)*a**4/n**3*x**3*mc2
    SE = pre*((4/3)*Lr*s - (4/3)*K(n, l) + (10/9)*s + (139/32 - 2*math.log(2))*math.pi*a*s)
    VP = pre*(-(4/15) + 5*math.pi/48*a)*s
    B40 = 1.5*math.pi**2*math.log(2) - 10*math.pi**2/27 - 2179/648 - 2.25*z3
    TL = (a/math.pi)**2*a**4/n**3*x**3*mc2*B40*s
    if l == 0: an = -2*(math.log(2/n) + sum(1/i for i in range(1, n + 1)) + 1 - 1/(2*n))
    else: an = 1/(l*(l + 1)*(2*l + 1))
    mN, me, mr = mpme, 1.0, x
    ES = x**3*r*a**5/(math.pi*n**3)*mc2*((1/3)*s*La - (8/3)*K(n, l) - (1/9)*s - (7/3)*an
                                         - 2*s*(mN**2*math.log(me/mr) - me**2*math.log(mN/mr))/(mN**2 - me**2))
    RR = x**3*r*a**6/(math.pi**2*n**3)*mc2*(6*z3 - 2*math.pi**2*math.log(2) + 35*math.pi**2/36 - 448/27 + (2/3)*math.pi*a*La**2)*s
    NS = (2/3)*x**3*a**2/n**3*mc2*(a*Rp/lamC)**2*s
    return dict(SE=SE, VP=VP, TL=TL, ES=ES, RR=RR, NS=NS)
def D(n, l): return sum(parts(n, l).values())
MHz_per_ueV = 241.79892623
out = ['hydrogen: alpha^-1 %.9f (engine), m_e c^2 %.6f MHz, m_p/m_e %.9f; x = m_r/m_e %.12f' % (1/a, mc2, mpme, x)]
for n, l in [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2), (4, 0), (4, 1)]:
    p = parts(n, l)
    out.append('import n %d l %d: ' % (n, l) + '  '.join('%s %+.4f' % (k, v) for k, v in p.items()) + '  total %+.4f MHz (%+.4f ueV)' % (D(n, l), D(n, l)/MHz_per_ueV))
# engine intervals before
F = {k: E[k] for k in E if k.startswith('F_')}
lamb0 = F['F_2S'] - F['F_2P1/2']; lamb = lamb0 + D(2, 0) - D(2, 1)
out.append('2S1/2 - 2P1/2 (spin centroid): engine %.4f -> predicted %.4f MHz vs 1057.8298(32) (Bezginov et al., Science 365, 1007 (2019)): %+.4f MHz;'
           ' vs 1057.845(9) (Lundeen, Pipkin, PRL 46, 232 (1981)): %+.4f MHz' % (lamb0, lamb, lamb - 1057.8298, lamb - 1057.845))
s3 = F['F_2P3/2'] - F['F_2S']; s3n = s3 + D(2, 1) - D(2, 0)
out.append('2P3/2 - 2S1/2: engine %.4f -> predicted %.4f MHz vs 9911.200(12) (Hagley, Pipkin, PRL 72, 1172 (1994)): %+.4f MHz' % (s3, s3n, s3n - 9911.200))
out.append('2P3/2 - 2P1/2: engine %.4f MHz, unchanged by the import (same Delta for both j)' % (F['F_2P3/2'] - F['F_2P1/2']))
c12 = 0.75*E['binding_MHz'] + F['F_2S'] - F['F_1S']; meas = 2466061413.187035
out.append('1S - 2S hyperfine centroid: engine %.4f MHz vs 2 466 061 413.187035(10) (Parthey et al., PRL 107, 203001 (2011), centroid): residual %+.4f -> %+.4f MHz'
           % (c12, c12 - meas, c12 - meas + D(2, 0) - D(1, 0)))
out.append('  (429 compared the F = 1 component of the model with this centroid: U(2S, F=1) - U(1S, F=1) = %+.4f MHz -> the quoted +6.8241 GHz; Parthey: centroid'
           ' = F = 1 line + 310.7122294 MHz)' % (E['U_2S_F1'] - E['U_1S_F1']))
Rinf_eng = a*a*mc2/2; Rinf_codata = 10973731.568157*299792458/1e6
out.append('  engine R_inf c %.4f MHz vs CODATA 2022 10 973 731.568 157(12) m^-1 x c = %.4f MHz: %+.4f MHz (%.2e), on 1S - 2S x (3/4) x m_r/m_e: %+.4f MHz'
           % (Rinf_eng, Rinf_codata, Rinf_eng - Rinf_codata, (Rinf_eng - Rinf_codata)/Rinf_codata, 0.75*x*(Rinf_eng - Rinf_codata)))
out.append('2p -> 1s line shift (both j): Delta(2P) - Delta(1S) = %+.4f MHz = %+.5f ueV' % (D(2, 1) - D(1, 0), (D(2, 1) - D(1, 0))/MHz_per_ueV))
out.append('3D - 3S: Delta(3D) - Delta(3S) = %+.4f MHz; 3D5/2 - 3S1/2 engine %.4f -> %.4f MHz vs 4013.204(78) (Clark et al. 1979, as in 516)'
           % (D(3, 2) - D(3, 0), F['F_3D5/2'] - F['F_3S'], F['F_3D5/2'] - F['F_3S'] + D(3, 2) - D(3, 0)))
if len(sys.argv) > 1:
    N = rd(sys.argv[1]); ok = True
    for k, (n, l) in {'F_1S': (1, 0), 'F_2S': (2, 0), 'F_2P1/2': (2, 1), 'F_2P3/2': (2, 1), 'F_3S': (3, 0), 'F_3P1/2': (3, 1), 'F_3D5/2': (3, 2),
                      'F_4S': (4, 0), 'F_4F7/2': (4, 3)}.items():
        d = N[k] - E[k]; ok &= abs(d - D(n, l)) < 0.01
        out.append('engine %s: new - old %+.4f vs import %+.4f MHz (diff %+.1e)' % (k, d, D(n, l), d - D(n, l)))
    for k in ('U_1S_F1', 'U_2S_F1', 'binding_MHz'): ok &= N[k] == E[k]
    out.append('H721a (engine level change = import within 0.01 MHz; U and binding unchanged): %s' % ok)
text = '\n'.join(out); print(text)
open(os.path.join(HERE, 'out', 'ana721.txt' if len(sys.argv) > 1 else 'pred721.txt'), 'w').write(text + '\n')

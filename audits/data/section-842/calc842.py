# Audit 842: even (S^2) part -- classical A^2/(2m) of the L_T form vs its spin-1/2 operator value and the Dirac FW hbar^2 term
# T_D = -hbar^2 (grad V)^2/(32 m^3 c^4) (Pachucki PRA 71, 012503, eq. 22).  Hydrogen test particle, metric of 808/838.  Written after prereg842.
import math, numpy as np, os, sys
D = os.path.dirname(os.path.abspath(__file__))
e0, me, hbar, c, eps0, h = 1.602176634e-19, 9.1093837139e-31, 1.054571817e-34, 299792458.0, 8.8541878188e-12, 6.62607015e-34
ge = 2.00231930436092; a_an = 0.5*(ge - 2.0); kappa = a_an + 0.5
ke = 1/(4*math.pi*eps0); k = ke*e0*e0; m = me; a0 = hbar**2/(m*k); a = 4*a0; n = math.sqrt(k/(m*a**3))
HEV = 4.135667696e-15
units = lambda f: f'{f/1e6:+.5f} MHz = {f*HEV*1e6:+.6f} ueV' + (f' (|lambda| {c/abs(f)*100:.4g} cm)' if f else '')
def gauss(N):
    x, w = np.polynomial.legendre.leggauss(N); u = 0.75 + 0.25*x[::-1]; return 1/u, 0.25*w[::-1]
NODES = dict(zip([f'u5{q}' for q in 'abcde'], zip(*gauss(5))))
nu_of = lambda Om: sum(w*L**2*Om[kk] for kk, (L, w) in NODES.items())/(2*math.pi)
LF = {'u5a': 0.5120091953982716, 'u5b': 0.5652161499007791, 'u5c': 0.66666666666666663, 'u5d': 0.8125025652578265, 'u5e': 0.9551918755393795}
def Efield(r): return r*(ke*e0/np.linalg.norm(r)**3)          # proton field at the electron
def omega(lf, dH, N):
    ecc = math.sqrt(max(0.0, 1 - lf*lf)); ec = math.sqrt(1 - ecc*ecc); P = np.array([1., 0, 0]); Q = np.array([0, 1., 0])
    s = 0.0; sw = 0.0
    for j in range(N):
        E = 2*math.pi*(j + 0.5)/N; w = 1 - ecc*math.cos(E)
        r = P*a*(math.cos(E) - ecc) + Q*a*ec*math.sin(E); v = (P*(-a*n*math.sin(E)) + Q*(a*n*ec*math.cos(E)))/w
        dx = 1e-6*np.linalg.norm(r); F = -np.array([(dH(r + d) - dH(r - d))/(2*dx) for d in np.eye(3)*dx])
        R = F/m; hv = np.cross(r, v); dA = np.cross(R, hv) + np.cross(v, np.cross(r, R))
        s += Q @ dA/(k*ecc/m)*w; sw += w
    return s/sw
CQ = kappa**2*hbar**2*k**2/(4*m**3*c**4); CT = -hbar**2*k**2/(32*m**3*c**4)
r4 = lambda lf: (1 + 0.5*(1 - lf*lf))/(a**4*lf**5)          # <r^-4>_t, 1 - e^2 = lf^2
out = {}
for N in (2048, 4096):
    for sk in ('H0', 'H1'):
        ln = [l.split() for l in open(f'{D}/out/params839.txt') if l.startswith(sk + ' ')][0]
        mu2 = np.array(list(map(float, ln[8:11]))); S = mu2/(ge*(-e0)/(2*m))
        dHc = lambda r, S=S: kappa**2*e0**2*np.dot(np.cross(Efield(r), S), np.cross(Efield(r), S))/(2*m**3*c**4)
        out[sk, 'cl', N] = nu_of({nd: omega(LF[nd], dHc, N) for nd in NODES})
    out['Q', N] = nu_of({nd: omega(LF[nd], lambda r: CQ/np.linalg.norm(r)**4, N) for nd in NODES})
    out['T', N] = nu_of({nd: omega(LF[nd], lambda r: CT/np.linalg.norm(r)**4, N) for nd in NODES})
anQ = CQ*(r4(1.0) - r4(0.5))/h; anT = CT*(r4(1.0) - r4(0.5))/h
print(f'stale: kappa {kappa:.8f}, a = 4 a0 {a:.6e} m, |S| = hbar/2; CQ/CT = {CQ/CT:.5f}')
okA = True
for sk, ref in (('H0', -0.0757e6), ('H1', -0.0758e6)):
    x = out[sk, 'cl', 4096]; okA &= abs(x - ref) <= 0.03*abs(ref)
    print(f'  {sk}: klasyczne A^2/(2m) nu {units(x)}  (841: {ref/1e6:+.4f}; N 2048 {out[sk, "cl", 2048]/1e6:+.5f})')
print(f"H842a (klasyczne A^2/2m = czesc parzysta 841 do 3 %): {'TAK' if okA else 'NIE'}")
okN = all(abs(out[q, 4096] - an) <= 0.01*abs(an) for q, an in (('Q', anQ), ('T', anT))) and \
      all(abs(out[kk, 4096] - out[kk, 2048]) <= 1e-3*abs(out[kk, 4096]) for kk in ('Q', 'T')) and \
      all(abs(out[sk, 'cl', 4096] - out[sk, 'cl', 2048]) <= 1e-3*abs(out[sk, 'cl', 4096]) for sk in ('H0', 'H1'))
print(f'  Q_Z\' (kwantowe A^2/2m): wezly {units(out["Q", 4096])}, analitycznie {units(anQ)}')
print(f'  T_D (Dirac hbar^2):      wezly {units(out["T", 4096])}, analitycznie {units(anT)}')
print(f"H842n (wezly = analityka do 1 %, N x 2 <= 0,1 %): {'TAK' if okN else 'NIE'}")
dm = out['Q', 4096] - out['T', 4096]
print(f"H842m (GLOWNA: |nu[Q_Z'] - nu[T_D]| <= 0,015 MHz): {'TAK' if abs(dm) <= 0.015e6 else 'NIE'}  (roznica {units(dm)})")
mc2a6 = m*c*c*(k/(hbar*c))**6/h
print(f'  poziom 2P (<r^-4> = 1/(24 a0^4)): Q_Z\' {CQ/(24*a0**4)/h/1e6:+.5f} MHz, T_D {CT/(24*a0**4)/h/1e6:+.5f} MHz '
      f'(= -m c^2 alpha^6 / 768 = {-mc2a6/768/1e6:+.5f})')
print(f'  wspolczynnik c przy q^2 hbar^2 E^2/(m^3 c^4): silnik kappa^2/4 = {kappa**2/4:+.6f}, Dirac -1/32 = {-1/32:+.6f} '
      f'(kappa_eff^2 = -1/8: brak rzeczywistego kappa)')

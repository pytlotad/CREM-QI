# Audit 438: exact resonant drive integral on a frozen Kepler orbit (hbar = mu = k = 1, theta = psi + J (E + e sin E))
# D(psi) per orbit = oint cos(theta) dg = -sin(psi) I + cos(psi) I_c, I = int g'(E) sin f dE, I_c = int g'(E) cos f dE
import math
def quad(fun,a,b,**kw):
    # composite Simpson, 2e5 intervals (integrand smooth and periodic -> spectrally accurate trapezoid would also do)
    n=200000; h=(b-a)/n; s=fun(a)+fun(b)
    for i in range(1,n): s+=(4 if i%2 else 2)*fun(a+i*h)
    return (s*h/3,0)
for J in (1.0, 2.0, 3.0, 2.5):
  for e in (0.3, 0.6, 0.9):
    for r0 in (2.0, 0.5):
      a = J*J
      g = lambda E: math.exp(-a*(1-e*math.cos(E))/r0)
      gp = lambda E: g(E)*(-a*e*math.sin(E)/r0)
      f = lambda E: J*(E+e*math.sin(E))
      I = quad(lambda E: gp(E)*math.sin(f(E)), -math.pi, math.pi, limit=400, epsabs=1e-15, epsrel=1e-13)[0]
      Ic = quad(lambda E: gp(E)*math.cos(f(E)), -math.pi, math.pi, limit=400, epsabs=1e-15, epsrel=1e-13)[0]
      print(f"J {J} e {e} r0 {r0}: I = {I:+.3e}  I_c = {Ic:+.3e}")

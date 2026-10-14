# Audit 438 design (exact quadrature, frozen Kepler orbit, hbar = mu = k = 1, theta = psi + J (E + e sin E)).
# Per unit eps, orbit-averaged:  dJ/dt = J^3 (1/T) oint cos(theta) dg = -(J^3/T) I sin(psi)        (I_c = 0 by parity)
#                                 clock pull term <(g - gbar) sin theta> = C sin(psi),  C = (1/T) oint (g - gbar) cos f dt
import math
def simpson(fun,a,b,n=20000):
    h=(b-a)/n; s=fun(a)+fun(b)
    for i in range(1,n): s+=(4 if i%2 else 2)*fun(a+i*h)
    return s*h/3
for r0 in (1.0,):
  for J in (2.0,3.0):
    e=0.6; a=J*J; T=2*math.pi*J**3
    g=lambda E: math.exp(-a*(1-e*math.cos(E))/r0)
    gp=lambda E: g(E)*(-a*e*math.sin(E)/r0)
    f=lambda E: J*(E+e*math.sin(E))
    dtdE=lambda E: (1-e*math.cos(E))*J**3          # dt/dE = r/(a w) = (1 - e cos E) J^3
    I=simpson(lambda E: gp(E)*math.sin(f(E)),-math.pi,math.pi)
    gbar=simpson(lambda E: g(E)*dtdE(E),-math.pi,math.pi)/T
    C=simpson(lambda E: (g(E)-gbar)*math.cos(f(E))*dtdE(E),-math.pi,math.pi)/T
    drive=J**3/T*I        # dJ/dt per eps = -drive sin psi
    w=1/J**3
    print(f"r0 {r0} J {J}: I {I:+.5e}  gbar {gbar:.5e}  C {C:+.5e}  drive amplitude dJ/dt/eps = {abs(drive):.5e}  omega {w:.5f}")

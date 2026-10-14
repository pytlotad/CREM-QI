# Audit 440: sign of C_k * I_k (stability of rung k for kappa > 0) on frozen Kepler orbits, hbar = mu = k = 1, J = k.
# I_k = int g'(E) sin[k(E + e sin E)] dE ; C_k = (1/T) oint (g - gbar) cos[k(E + e sin E)] dt ; stable for kappa > 0 iff C_k I_k > 0
# (438: k = 2, r0 = 1: I = -0,130, C = -0,0108 -> product > 0, stable; k = 3: I = -0,017, C = +5,3e-4 -> unstable).
import math
def simpson(fun,a,b,n=8000):
    h=(b-a)/n; s=fun(a)+fun(b)
    for i in range(1,n): s+=(4 if i%2 else 2)*fun(a+i*h)
    return s*h/3
shapes={'exp r0=0.5':lambda r:math.exp(-r/0.5),'exp r0=1':lambda r:math.exp(-r/1.0),'exp r0=2':lambda r:math.exp(-r/2.0),
        'exp r0=4':lambda r:math.exp(-r/4.0),'1/r':lambda r:1/r,'1/r^2':lambda r:1/r**2,'1/r^3':lambda r:1/r**3}
for name,g in shapes.items():
  line=f"{name:11s}"
  for e in (0.3,0.6,0.9):
    signs=""
    for k in range(1,7):
      J=float(k); a=J*J; T=2*math.pi*J**3
      r=lambda E: a*(1-e*math.cos(E)); dtdE=lambda E:(1-e*math.cos(E))*J**3; f=lambda E: J*(E+e*math.sin(E))
      h=1e-6*a; gp=lambda E: (g(r(E)+h)-g(r(E)-h))/(2*h)*a*e*math.sin(E)
      I=simpson(lambda E: gp(E)*math.sin(f(E)),-math.pi,math.pi)
      gbar=simpson(lambda E: g(r(E))*dtdE(E),-math.pi,math.pi)/T
      C=simpson(lambda E:(g(r(E))-gbar)*math.cos(f(E))*dtdE(E),-math.pi,math.pi)/T
      p=C*I; signs+=("+" if p>0 else ("-" if p<0 else "0"))
    line+=f"  e={e}: {signs}"
  print(line)

# POST HOC 817: split nu_815 - nu_true into dynamics (engine vs Darwin Lagrangian at the start invariants) and start label.
# Usage (from section-817): python3 posthoc/split817.py > out/posthoc817.txt
import re, math, numpy as np
x,w=np.polynomial.legendre.leggauss(5); u=0.75+0.25*x[::-1]; N=dict(zip([f'u5{c}' for c in 'abcde'],zip(1/u,0.25*w[::-1])))
Oe,Od={},{}
for ln in open('out/ana817.txt'):
    m=re.match(r'\s+(u5\w): Omega_silnik\(eps -> 0\) (\S+), Omega_Darwin\(start\) (\S+)',ln)
    if m and m.group(1) not in Oe: Oe[m.group(1)]=float(m.group(2)); Od[m.group(1)]=float(m.group(3))
nu=lambda O: sum(ww*L*L*O[k] for k,(L,ww) in N.items())/(2*math.pi)
nt=5474.7649e6
for k in N: print(f'{k}: Omega_silnik - Omega_Darwin = {Oe[k]-Od[k]:.4e} rad/s; wzgl. {Oe[k]/Od[k]-1:+.3e} = {(Oe[k]/Od[k]-1)*137.035999177**2:+.4f} alfa^2')
print(f'nu(Omega_silnik) = {nu(Oe)/1e6:.4f} MHz; nu(Omega_Darwin(start)) = {nu(Od)/1e6:.4f} MHz; nu_true = 5474.7649 MHz')
print(f'czesc dynamiki (silnik - Darwin na startach) {(nu(Oe)-nu(Od))/1e6:+.4f} MHz; czesc etykiety (Darwin na startach - torusy) {(nu(Od)-nt)/1e6:+.4f} MHz')

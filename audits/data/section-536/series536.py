# Audit 536c (sympy): QM P_S = y^2/(1+y^2), y = x/(1+sqrt(1+x^2)), vs classical stationary w = (x/2)^2 at K = dE/4.
import sympy as sp
x=sp.symbols('x',positive=True)
y=x/(1+sp.sqrt(1+x**2)); PS=y**2/(1+y**2); closed=(1-1/sp.sqrt(1+x**2))/2
print("P_S == (1 - 1/sqrt(1+x^2))/2:",sp.simplify(PS-closed)==0)
print("P_S series:",sp.series(PS,x,0,8))
cl=(x/2)**2
print("classical w - P_S:",sp.series(cl-PS,x,0,8))
print("relative excess (w - P_S)/P_S:",sp.series((cl-PS)/PS,x,0,6))

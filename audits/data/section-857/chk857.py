# Audit 857: check of the closed form against direct integration of da/dt, de/dt (units C = 1, a0 = 1; no engine).
import math
def closed(a0, e0, r):
    piece = lambda a, e: a**3 * (1 - e*e)**2.5 / (1 + math.sqrt(1 - e*e))**2
    ef = 0.0
    if e0 > 0:
        A = a0 * (1 - e0*e0) / e0**(4/3); lo, hi = 0.0, e0
        for _ in range(200):
            m = 0.5*(lo + hi); (lo, hi) = (m, hi) if A*m**(4/3)/(1 - m*m) < r else (lo, m)
        ef = 0.5*(lo + hi)
    return 4/3*(piece(a0, e0) - piece(r, ef))
def ode(a0, e0, r, N=400000):
    # integrate in u = ln a from ln a0 down to ln r with RK4 on (t, e)
    def f(u, y):
        a = math.exp(u); t, e = y; dadt = -(1 + e*e/2)/(a*a*(1 - e*e)**2.5); dedt = -0.75*e/(a**3*(1 - e*e)**1.5)
        return (a/dadt, a*dedt/dadt)
    u, y = math.log(a0), (0.0, e0); h = (math.log(r) - u)/N
    for _ in range(N):
        k1 = f(u, y); k2 = f(u + h/2, tuple(y[i] + h/2*k1[i] for i in range(2))); k3 = f(u + h/2, tuple(y[i] + h/2*k2[i] for i in range(2))); k4 = f(u + h, tuple(y[i] + h*k3[i] for i in range(2)))
        y = tuple(y[i] + h/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(2)); u += h
    return y[0]
for e0, r in ((0.0, 0.005), (0.3, 0.005), (0.6614378277661477, 0.00125), (0.8660254037844386, 0.005)):
    c, o, circ = closed(1, e0, r), ode(1, e0, r), (1 - r**3)/3
    print(f'CHK857 e0 {e0:.6f} r {r} closed {c:.10e} ode {o:.10e} rel {c/o - 1:+.2e}  closed/circular {c/circ:.6f}')

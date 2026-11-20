# Audit 815: independent integration of the engine's conservative e+e- equations WITH the Plummer separation floor (written after
# prereg815, before the runs).  Transcribed from pair_geometry.hpp / electrodynamics.hpp:
#   clampedSeparationVector: d_c = d * r_c/r, r_c = sqrt(r^2 + eps^2), eps = scale * r*,  r* = g hbar/(4 m c)
#   coulombForces: F_c = d_c * (kq/r_c^3) carried back by J^T (radial x r/r_c)  ->  F_1 = d * kq / r_c^3   (Plummer gradient)
#   darwinForces: darwinForceOnFirst(v1, v2, F_2/m2, d_c, kq)  (softened vector used directly; leading accelerations = softened Coulomb)
# Rest as pred814 (relativistic momenta, DOP853, apocentres from (x1 - x2).(v1 - v2) = 0, Omega = slope of apse angle on time).
import math, sys, os
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-814'))
from pred814 import parse_start, darwin_on_first
HBAR = 6.62607015e-34/(2*math.pi); G_E = 2.00231930436092; ME = 9.1093837139e-31; CL = 299792458.0
RSTAR = G_E*HBAR/(4*ME*CL)

def omega(st, scale, orbits=30, rtol=1e-13):
    m1, m2, kq, c = st['m1'], st['m2'], st['kq'], st['c']; eps = scale*RSTAR
    L0 = np.linalg.norm(st['x1'] - st['x2']); T0 = st['period']; P0 = m1*np.linalg.norm(st['v1'] - st['v2'])
    vel = lambda p, m: p/(m*math.sqrt(1 + (p @ p)/(m*m*c*c)))
    mom = lambda v, m: v*m/math.sqrt(1 - (v @ v)/(c*c))
    def rhs(t, y):
        x1, x2, p1, p2 = y[0:3]*L0, y[3:6]*L0, y[6:9]*P0, y[9:12]*P0
        v1, v2 = vel(p1, m1), vel(p2, m2); d = x1 - x2; r = np.linalg.norm(d); rc = math.sqrt(r*r + eps*eps)
        dc = d*(rc/r)
        F1 = d*(kq/rc**3); F2 = -F1
        D1 = darwin_on_first(v1, v2, F2/m2, dc, kq, c); D2 = darwin_on_first(v2, v1, F1/m1, -dc, kq, c)
        return np.concatenate([v1*T0/L0, v2*T0/L0, (F1 + D1)*T0/P0, (F2 + D2)*T0/P0])
    def ev(t, y):
        x1, x2, p1, p2 = y[0:3]*L0, y[3:6]*L0, y[6:9]*P0, y[9:12]*P0
        return (x1 - x2) @ (vel(p1, m1) - vel(p2, m2))/(L0*L0/T0)
    ev.direction = -1
    y0 = np.concatenate([st['x1']/L0, st['x2']/L0, mom(st['v1'], m1)/P0, mom(st['v2'], m2)/P0])
    sol = solve_ivp(rhs, (0, orbits + 0.5), y0, method='DOP853', rtol=rtol, atol=rtol*1e-3, events=ev)
    te = sol.t_events[0]; ye = sol.y_events[0]
    d0 = st['x1'] - st['x2']; vr0 = st['v1'] - st['v2']; Lh = np.cross(d0, vr0); Lh /= np.linalg.norm(Lh)
    e1 = d0/np.linalg.norm(d0); e2 = np.cross(Lh, e1)
    ph = np.unwrap([math.atan2((y[0:3] - y[3:6]) @ e2, (y[0:3] - y[3:6]) @ e1) for y in ye])
    t = te*T0
    return np.polyfit(t, ph, 1)[0], len(t)

if __name__ == '__main__':   # smoke: pred(0) must equal pred814 on the 814 start
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-814'))
    import pred814
    st = parse_start(sys.argv[1])
    for s in [float(x) for x in sys.argv[2:]]:
        print(s, '%.10e %d' % omega(st, s))
    print('pred814', '%.10e' % pred814.omega(st)[0])

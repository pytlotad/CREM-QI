# Audit 814: INDEPENDENT integration of the engine's conservative pair equations (written before the runs).
#   p_i = gamma m v_i;  dp_1/dt = F_C + F_D(v1, v2, a2_lead, x1 - x2),  F_C = kq1q2 (x1 - x2)/r^3,  a_lead = Coulomb force / m (non-relativistic)
#   F_D = Euler-Lagrange force of L_D = C [v1.v2 + (v1.n)(v2.n)], C = kq1q2/(2 c^2 r), accelerations reduced to Coulomb (electrodynamics.hpp
#   darwinForceOnFirst, transcribed term by term).  Start = the engine's PAIR CONS814 line (an input, not a result).
#   Apocentres: event (x1 - x2).(v1 - v2) = 0 crossing + -> -; apse angle in the start orbital plane; Omega = linear regression on time.
import math
import numpy as np
from scipy.integrate import solve_ivp

def parse_start(path, tag='mechB0'):
    for ln in open(path):
        if ln.startswith(f'PAIR CONS814 {tag} start'):
            w = ln.split(); g = lambda k: float(w[w.index(k) + 1])
            i = w.index('period') + 2; v = [float(x) for x in w[i:i + 12]]
            return dict(m1=g('m1'), m2=g('m2'), K=g('K'), kq=g('kq1q2'), c=g('c'), dt=g('dt'), period=g('period'),
                        x1=np.array(v[0:3]), x2=np.array(v[3:6]), v1=np.array(v[6:9]), v2=np.array(v[9:12]))
    return None

def darwin_on_first(v1, v2, a2, d, kq, c):
    r = np.linalg.norm(d); n = d/r
    r1 = v1 @ n; r2 = v2 @ n; rr = r1 - r2
    rel = v1 - v2; nrate = (rel - n*rr)/r
    r2rate = a2 @ n + v2 @ nrate
    C = kq/(2*c*c*r); vp = v1 @ v2
    spatial = (v2*r1 + v1*r2 - n*(vp + 3*r1*r2))*(C/r)
    mom = (v2 + n*r2)*(-C*rr/r) + (a2 + nrate*r2 + n*r2rate)*C
    return spatial - mom

def omega(st, orbits=30, rtol=1e-13):
    m1, m2, kq, c = st['m1'], st['m2'], st['kq'], st['c']
    L0 = np.linalg.norm(st['x1'] - st['x2']); T0 = st['period']; P0 = m1*np.linalg.norm(st['v1'] - st['v2'])
    vel = lambda p, m: p/(m*math.sqrt(1 + (p @ p)/(m*m*c*c)))
    mom = lambda v, m: v*m/math.sqrt(1 - (v @ v)/(c*c))
    def rhs(t, y):
        x1, x2, p1, p2 = y[0:3]*L0, y[3:6]*L0, y[6:9]*P0, y[9:12]*P0
        v1, v2 = vel(p1, m1), vel(p2, m2); d = x1 - x2; r = np.linalg.norm(d)
        F1 = d*(kq/r**3); F2 = -F1
        D1 = darwin_on_first(v1, v2, F2/m2, d, kq, c); D2 = darwin_on_first(v2, v1, F1/m1, -d, kq, c)
        return np.concatenate([v1*T0/L0, v2*T0/L0, (F1 + D1)*T0/P0, (F2 + D2)*T0/P0])
    def ev(t, y):
        x1, x2, p1, p2 = y[0:3]*L0, y[3:6]*L0, y[6:9]*P0, y[9:12]*P0
        return (x1 - x2) @ (vel(p1, m1) - vel(p2, m2))/(L0*L0/T0)
    ev.direction = -1
    y0 = np.concatenate([st['x1']/L0, st['x2']/L0, mom(st['v1'], m1)/P0, mom(st['v2'], m2)/P0])
    sol = solve_ivp(rhs, (0, orbits + 0.5), y0, method='DOP853', rtol=rtol, atol=rtol*1e-3, events=ev, dense_output=False)
    te = sol.t_events[0]; ye = sol.y_events[0]
    d0 = st['x1'] - st['x2']; vr0 = st['v1'] - st['v2']; Lh = np.cross(d0, vr0); Lh /= np.linalg.norm(Lh)
    e1 = d0/np.linalg.norm(d0); e2 = np.cross(Lh, e1)
    ph = np.unwrap([math.atan2((y[0:3] - y[3:6]) @ e2, (y[0:3] - y[3:6]) @ e1) for y in ye])
    t = te*T0; k = np.polyfit(t, ph, 1)[0]   # successive apocentres differ by the small advance delta (mod 2 pi): slope = Omega_aps
    Tr = np.polyfit(np.arange(len(t)), t, 1)[0]
    return k, len(t), Tr

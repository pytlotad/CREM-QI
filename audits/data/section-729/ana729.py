"""Audit 729 (offline): spin length S + 1/2 from the single-valuedness of the geometric phase (Dirac) and Maslov on the cap; does the 728 phase rule
carry over to spin-spin coupling (Ps singlet / triplet)?  See prereg729.txt."""
import math, os
import numpy as np
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__))
out = []
def P(s): out.append(s); print(s, flush=True)
def wrap(x): return (x + math.pi) % (2*math.pi) - math.pi


def berry(n, S):
    """n: unit vectors along a closed loop (3, M); north / south gauge geometric phase S (cos th -+ 1) dphi, trapezoid"""
    th = np.arccos(np.clip(n[2], -1, 1)); ph = np.unwrap(np.arctan2(n[1], n[0]))
    c = np.cos(th); dph = np.diff(ph); cm = 0.5*(c[1:] + c[:-1])
    return S*np.sum((cm - 1)*dph), S*np.sum((cm + 1)*dph)


# ---- N1
P('N1: loop around an axis tilted 0.3 rad from z, cone 0.5 rad (encircles the north pole); geometric phase in north / south gauge')
a = np.array([math.sin(0.3), 0, math.cos(0.3)]); e1 = np.array([math.cos(0.3), 0, -math.sin(0.3)]); e2 = np.array([0, 1.0, 0])
phi = np.linspace(0, 2*np.pi, 400001)
n = (math.cos(0.5)*a)[:, None] + math.sin(0.5)*(np.cos(phi)*e1[:, None] + np.sin(phi)*e2[:, None])
okA = True
for lab, S in (('1/2', 0.5), ('1', 1.0), ('3/2', 1.5), ('2', 2.0), ('sqrt 2', math.sqrt(2)), ('sqrt(3)/2', math.sqrt(3)/2)):
    bn, bs = berry(n, S); diff = bs - bn; single = abs(wrap(diff)) < 1e-8
    okA &= abs(diff/(4*np.pi*S) - 1) < 1e-8 and (single == (abs(2*S - round(2*S)) < 1e-12))
    P('  |S| = %-9s: north %+.10f, south %+.10f, difference / (4 pi |S|) = %.12f, wrapped difference %+.2e -> single-valued: %s' %
      (lab, bn, bs, diff/(4*np.pi*S), wrap(diff), single))
P('H729a (difference = 4 pi |S|; single-valued iff 2|S| integer): %s' % okA)

# ---- N2
P('\nN2: precession loop about z at polar angle theta; turning points of S_x; cap quantization with / without Maslov')
th0 = 1.1; phi = np.linspace(0.1, 0.1 + 2*np.pi*3, 600001); sx = np.sin(th0)*np.cos(phi); d = np.diff(sx)   # start off an extremum (run 1 started at phi = 0, an S_x maximum)
tp = int(np.sum(np.sign(d[1:]) != np.sign(d[:-1])))/3
P("  turning points of S_x per loop (3 loops sampled): %.3f" % tp)
okB = tp == 2
for lab, S, want in (('1/2', 0.5, [0.0]), ('1', 1.0, [0.5, -0.5]), ('3/2', 1.5, [1.0, 0.0, -1.0]), ('2', 2.0, [1.5, 0.5, -0.5, -1.5]), ('sqrt 2', math.sqrt(2), None)):
    withM = [S - 0.5 - k for k in range(0, 20) if abs(S - 0.5 - k) < S]
    without = [S - k for k in range(1, 20) if abs(S - k) < S]
    if want is not None: okB &= np.allclose(withM, want)
    P('  |S| = %-6s: with Maslov S_z = %s (%d states);  without: %s (%d)' % (lab, ', '.join('%+.3f' % x for x in withM), len(withM),
                                                                         ', '.join('%+.3f' % x for x in without), len(without)))
P('H729b (2 turning points; with Maslov: |S| = 1 -> +-1/2, 3/2 -> 0, +-1, 1/2 -> 0; 2|S| states): %s' % okB)

# ---- N3
P('\nN3: two spins |s1| = |s2| = 1 coupled by A s1.s2 (A = 0.01): geometric phase per precession period about S^ (north gauge), DOP853')
A = 0.01; okC = True
for Sl in (0.5, 1.0, 1.5, 2.0 - 1e-9):
    ca = Sl/2; sa = math.sqrt(max(0.0, 1 - ca*ca))
    s1 = np.array([sa, 0, ca]); s2 = np.array([-sa, 0, ca])
    rot = np.array([[1, 0, 0], [0, math.cos(0.4), -math.sin(0.4)], [0, math.sin(0.4), math.cos(0.4)]])   # generic orientation of S
    s1, s2 = rot @ s1, rot @ s2
    def f(t, y): x1, x2 = y[:3], y[3:]; return np.concatenate([A*np.cross(x2, x1), A*np.cross(x1, x2)])
    T = 2*np.pi/(A*Sl)
    sol = solve_ivp(f, (0, T), np.concatenate([s1, s2]), method='DOP853', rtol=1e-12, atol=1e-14, dense_output=True)
    ts = np.linspace(0, T, 200001); Y = sol.sol(ts)
    Sh = (s1 + s2)/np.linalg.norm(s1 + s2)
    ez = Sh; ex = np.cross([0, 1.0, 0], ez); ex /= np.linalg.norm(ex); ey = np.cross(ez, ex)
    tot = 0.0; parts = []
    for k in (0, 3):
        v = Y[k:k + 3]; nv = np.stack([(v*ex[:, None]).sum(0), (v*ey[:, None]).sum(0), (v*ez[:, None]).sum(0)])
        bn, _ = berry(nv/np.linalg.norm(nv, axis=0), 1.0); tot += bn; parts.append(bn)
    pred = 2*np.pi*(Sl - 2)
    okC &= abs(wrap(tot - pred)) < 1e-7
    P('  |S| = %.3f: phases %+.9f, %+.9f; sum %+.9f vs 2 pi (|S| - 2) %+.9f (diff %+.1e); sum mod 2 pi = %+.6f -> single-valued: %s' %
      (Sl, parts[0], parts[1], tot, pred, wrap(tot - pred), wrap(tot), abs(wrap(tot)) < 1e-6))
P('H729c (sum = 2 pi (|S| - 2); 728 rule selects |S| = 1, 2, rejects rule-Q 1/2, 3/2): %s' % okC)

# ---- H729d analytic table
P('\nH729d: triangle-parity rule |J1| + |J2| + |J| - 3/2 integer (phase 2 pi (|J| - |J1| - |J2|) - pi per precession cycle) vs the 728 rule')
rows = [('o-Ps orbit + spin', 1.5, 1.5, [0.5, 1.0, 1.5, 2.0, 2.5]), ('hydrogen orbit + spin', 1.5, 1.0, [0.5, 1.0, 1.5, 2.0, 2.5]),
        ('Ps spin + spin', 1.0, 1.0, [0.5, 1.0, 1.5, 2.0])]
okD = True
for lab, J1, J2, cands in rows:
    tri = [J for J in cands if abs(J1 - J2) <= J <= J1 + J2]
    par = [J for J in tri if abs((J1 + J2 + J - 1.5) - round(J1 + J2 + J - 1.5)) < 1e-12]
    r728 = [J for J in tri if abs((J - J2) - round(J - J2)) < 1e-12] if 'orbit' in lab else [J for J in tri if abs((J - J1 - J2) - round(J - J1 - J2)) < 1e-12]
    ruleQ = {'o-Ps orbit + spin': [0.5, 1.5, 2.5], 'hydrogen orbit + spin': [1.0, 2.0], 'Ps spin + spin': [0.5, 1.5]}[lab]
    okD &= par == ruleQ
    P('  %-22s |J1| %.1f |J2| %.1f: parity rule %s; 728 rule %s; engine rule Q %s' % (lab, J1, J2, par, r728, ruleQ))
P('H729d (parity rule = rule Q in all three; 728 rule = rule Q for orbit + spin only): %s' % okD)
open(os.path.join(HERE, 'out', 'ana729.txt'), 'w').write('\n'.join(out) + '\n')

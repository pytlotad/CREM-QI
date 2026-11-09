# Audit 762 (offline): <r^-3> from the EBK wave at order hbar^2 (737 iteration) -> fine structure from spin precession.  -> out/ana762.txt
import os, sys, io, contextlib, math, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(D, '..', 'section-737'))
with contextlib.redirect_stdout(io.StringIO()): import ana737 as w
out = []; pr = lambda s: (print(s, flush=True), out.append(s))
ALPHA = 7.2973525643e-3; AE = 1.15965218059e-3; MP_ME = 1836.15267343; MU = 1 / (1 + 1 / MP_ME)
HZ = 6.579683920502e15; EV = 27.211386245988; c_cm = 2.99792458e10; MEAS = 9911.200 + 1057.845; XI_PS = 2741.682
un = lambda f: f"{f:.3f} MHz = {f*1e6*EV/HZ*1e6:.4f} µeV = {c_cm/(f*1e6):.4f} cm"
res = {}
for n, l in [(2, 1), (3, 1), (3, 2), (4, 1), (4, 2), (4, 3)]:
    r = w.grid(n); key = (n, 1); q = 1 / (n**3 * l * (l + 0.5) * (l + 1)); uq = w.u_qm(n, l, r); r3q = np.trapezoid(uq * uq / r**3, r)
    v = {}
    for k, g in [(0, False), (1, False), (1, True)]:
        u, ph = w.ebk(n, l, r, k, key, g); v[(k, g)] = np.trapezoid(u * u / r**3, r)
    w._cache.clear(); res[(n, l)] = (q, v)
    pr(f"{n}{'spdf'[l]}: QM {q:.8f} (siatka {r3q/q-1:+.1e}); k0 {v[(0,False)]/q-1:+.3e}; k1-in {v[(1,False)]/q-1:+.3e}; k1-gl {v[(1,True)]/q-1:+.3e}")
q, v = res[(2, 1)]
def fs(r3, g=2 * (1 + AE)): return 0.5 * ALPHA**2 * r3 * 1.5 * MU**3 * ((g - 1) + g / MP_ME) * HZ / 1e6
pr(f"pomiar 2P3/2 - 2P1/2: {un(MEAS)}")
for name, kk in [('k0 (733)', (0, False)), ('k1-in', (1, False)), ('k1-gl', (1, True))]:
    f = fs(v[kk]); pr(f"  {name}: {un(f)}; - pomiar {f-MEAS:+.3f} MHz ({(f/MEAS-1)*100:+.4f} %)")
    xi = (3 / 16) * ALPHA**2 * v[kk] * ((2 * 2 * (1 + AE) - 1) / 3) * HZ / 1e6
    pr(f"    o-Ps 2^3P_J stala LS: {xi:.3f} MHz wobec {XI_PS} z teorii QED ({(xi/XI_PS-1)*100:+.4f} %); p-Ps 2^1P_1: bez LS")
pr(f"  QM <r^-3>: {un(fs(q))}")
a = abs(v[(1, False)] / q - 1) < 1e-3 and abs(fs(v[(1, False)]) / MEAS - 1) < 1e-3
b = all(abs(vv[(1, False)] / qq - 1) < abs(vv[(0, False)] / qq - 1) for qq, vv in res.values())
pr(f"H762a: {'TAK' if a else 'NIE'};  H762b: {'TAK' if b else 'NIE'}")
open(os.path.join(D, 'out', 'ana762.txt'), 'w').write('\n'.join(out) + '\n')

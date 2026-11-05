# Audit 717 (offline): the relativistic joint clock resonance 2 D = j_r hbar Omega_r + j_phi hbar Omega_phi (715) at order alpha^4, against
# Sommerfeld / Dirac / EBK and the measured fine structure (hydrogen) and n = 2 intervals (positronium, per channel).  Atomic units, c = 1/alpha.
# usage: python3 ana717.py
import os, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ALPHA = 7.2973525693e-3; C = 1/ALPHA; C2 = C*C
HA_HZ = 6.579683920502e15; HZ_EV = 4.135667696e-15; CL = 299792458.0
A4 = ALPHA**2*HA_HZ           # alpha^4 m c^2 in Hz
out = []


def units(hz): return f'{hz/1e6:+,.2f} MHz = {hz*HZ_EV*1e6:+.4f} ueV ({CL/abs(hz)*100:.4g} cm)' if hz else '0'


class System:
    def __init__(self, ps): self.ps = ps; self.mu = 0.5 if ps else 1.0

    def H(self, r, p, pr):
        if self.ps: return 2*(C2*math.sqrt(1 + p*p/C2) - C2) - 1/r - (p*p + pr*pr)/(2*C2*r)
        return C2*math.sqrt(1 + p*p/C2) - C2 - 1/r

    def rhs(self, t, y):
        x, yy, px, py = y[:4]; r2 = x*x + yy*yy; r = math.sqrt(r2); p2 = px*px + py*py; pr = (x*px + yy*py)/r
        g = math.sqrt(1 + p2/C2)
        if self.ps:
            vx = 2*px/g - (px + pr*x/r)/(C2*r); vy = 2*py/g - (py + pr*yy/r)/(C2*r)
            pdr = x*px + yy*py; Q = p2*r2 + pdr*pdr
            fx = -x/r**3 + (2*p2*x + 2*pdr*px)/(2*C2*r**3) - 3*Q*x/(2*C2*r**5)
            fy = -yy/r**3 + (2*p2*yy + 2*pdr*py)/(2*C2*r**3) - 3*Q*yy/(2*C2*r**5)
            u = (vx*vx + vy*vy)/4/C2; D = 2*C2*u/(1 + math.sqrt(1 - u))
        else:
            vx = px/g; vy = py/g; fx = -x/r**3; fy = -yy/r**3
            u = (vx*vx + vy*vy)/C2; D = C2*u/(1 + math.sqrt(1 - u))
        return [vx, vy, fx, fy, D, px*vx + py*vy]

    def period(self, E, L):
        a = 1/(2*abs(E)); e = math.sqrt(max(1e-12, 1 - L*L/(self.mu*a))); rp0 = a*(1 - e); ra0 = a*(1 + e)
        rp = brentq(lambda r: self.H(r, L/r, 0.0) - E, 0.5*rp0, math.sqrt(rp0*ra0), xtol=1e-15, rtol=1e-15)
        y0 = [rp, 0.0, 0.0, L/rp, 0.0, 0.0]
        T0 = 2*math.pi*math.sqrt(self.mu*a**3)
        tol = dict(method='DOP853', rtol=1e-13, atol=[1e-14*a, 1e-14*a, 1e-14, 1e-14, 1e-16, 1e-16])
        s1 = solve_ivp(self.rhs, (0, 0.5*T0), y0, **tol)
        ev = lambda t, y: y[0]*y[2] + y[1]*y[3]; ev.direction = 1; ev.terminal = True
        s2 = solve_ivp(self.rhs, (0.5*T0, 1.5*T0), s1.y[:, -1], events=ev, **tol)
        T = s2.t_events[0][0]; yT = s2.y_events[0][0]
        dphi = 2*math.pi + math.atan2(yT[1], yT[0])
        return dict(T=T, Wr=2*math.pi/T, Wphi=dphi/T, D=yT[4]/T, Jr=(yT[5] - L*dphi)/(2*math.pi), dphi=dphi)


def solve(sys_, L, cond, E0):
    return brentq(lambda E: cond(sys_.period(E, L)), E0*(1 + 2e-3), E0*(1 - 2e-3), xtol=1e-17, rtol=1e-15)


def part_H():
    S = System(False); ok = True
    out.append('PART 1: hydrogen (infinite nucleus), numerical resonance vs analytic (energies in alpha^4 m c^2 above Bohr)')
    def R0(n):
        N = brentq(lambda N: 2*C2*(1 - (1 + ALPHA**2/N**2)**-0.5) - n*(1/N**3)*(1 + ALPHA**2/N**2)**-1.5, n - 0.01, n + 0.01, xtol=1e-16)
        return C2*((1 + ALPHA**2/N**2)**-0.5 - 1)
    def som(n, L):
        N = (n - L) + math.sqrt(L*L - ALPHA**2); return C2*((1 + ALPHA**2/N**2)**-0.5 - 1)
    a4 = lambda E, n: (E + 0.5/n**2)/ALPHA**2
    for n in (1, 2):
        ref = R0(n); res = []
        out.append(f'  R0 analytic n {n}: {a4(ref, n):+.8f} (formula -3/8/n^4 = {-0.375/n**4:+.8f})')
        for L in (0.5, 1.0, 1.5):
            if L >= n: continue
            E = solve(S, L, lambda q: 2*q['D'] - n*q['Wr'], -0.5/n**2); res.append(E)
            eb = solve(S, L, lambda q: q['Jr'] - (n - L), -0.5/n**2)
            dev = (E - ref)/(ALPHA**2*0.375/n**4); ok &= abs(dev) < 1e-3
            out.append(f'    L {L}: E_res {a4(E, n):+.8f} (dev {dev:+.1e} of the alpha^4 term) | EBK J_r = n - L {a4(eb, n):+.8f} vs Sommerfeld {a4(som(n, L), n):+.8f}')
        if len(res) > 1: ok &= (max(res) - min(res))/(ALPHA**2*0.375/n**4) < 1e-3
    # R1: n = 2, k = 1, L = 1: j_r = 3, j_phi = -1
    E1 = solve(S, 1.0, lambda q: 2*q['D'] - (3*q['Wr'] - 1*q['Wphi']), -0.125)
    pred = som(2, 1.0) - 0.75*ALPHA**2/16
    dev = (E1 - pred)/(ALPHA**2*0.75/16); ok &= abs(dev) < 1e-3
    out.append(f'  R1 (n 2, k 1, L 1, j_r 3, j_phi -1): {a4(E1, 2):+.8f} vs Sommerfeld(2, 1) - 3/4/n^4 {a4(pred, 2):+.8f} (dev {dev:+.1e})')
    out.append(f'H717a (numerical = analytic within 1e-3 of the alpha^4 term; R0 independent of L): {ok}')


def part_fs():
    out.append('PART 2: hydrogen n = 2 and 1S - 2S at order alpha^4 (infinite-nucleus coefficients; Bohr 1S - 2S with reduced mass)')
    dirac = lambda n, j: -(n/(j + 0.5) - 0.75)/(2*n**4)
    dar = lambda n: 1/(2*n**3)
    so = lambda n, l, j: (j*(j + 1) - l*(l + 1) - 0.75)/(n**3*l*(l + 0.5)*(l + 1))/4 if l else 0.0   # (E_n^2/mc^2) n [...]/(l(l+1/2)(l+1)) / alpha^4
    R0 = lambda n: -0.375/n**4
    lv = {}
    lv['Dirac'] = {'1S': dirac(1, .5), '2S': dirac(2, .5), '2P1/2': dirac(2, .5), '2P3/2': dirac(2, 1.5)}
    lv['R0 A (spin via k)'] = {k: R0(int(k[0])) for k in lv['Dirac']}
    lv['R0 B (+ Darwin + SO)'] = {'1S': R0(1) + dar(1), '2S': R0(2) + dar(2), '2P1/2': R0(2) + so(2, 1, .5), '2P3/2': R0(2) + so(2, 1, 1.5)}
    lv['R1 (Dirac in n - 3/4/n^4)'] = {k: lv['Dirac'][k] - 0.75/int(k[0])**4 for k in lv['Dirac']}
    bohr = 2466038423.687e6; meas = {'Lamb': 1057.845e6, 'FS': 10969.045e6, '1S-2S': 2466061413.187e6}
    res = {}
    for name, e in lv.items():
        lamb = (e['2S'] - e['2P1/2'])*A4; fs = (e['2P3/2'] - e['2P1/2'])*A4; t12 = bohr + (e['2S'] - e['1S'])*A4
        res[name] = (lamb, fs, t12)
        out.append(f'  {name:26s}: 2S1/2 - 2P1/2 {units(lamb)} | 2P3/2 - 2P1/2 {units(fs)} | 1S - 2S {t12/1e6:,.1f} MHz (- measured {(t12 - meas["1S-2S"])/1e6:+,.1f})')
    out.append(f'  measured: 2S1/2 - 2P1/2 1,057.845 MHz (Lundeen-Pipkin 1981), 2P3/2 - 2P1/2 10,969.045 MHz (Hagley-Pipkin 1994 + LP), 1S - 2S 2,466,061,413.187 MHz (Parthey 2011)')
    A = res['R0 A (spin via k)']; B = res['R0 B (+ Darwin + SO)']; R1 = res['R1 (Dirac in n - 3/4/n^4)']
    out.append(f'H717b (R0, reading A: 2P3/2 - 2P1/2 = 0): {abs(A[1]) < 1.0}')
    out.append(f'H717c (R0 B: 2S1/2 - 2P1/2 > 20x measured; R1: |1S - 2S - measured| > 100 GHz): {B[0] > 20*meas["Lamb"] and abs(R1[2] - meas["1S-2S"]) > 100e9}'
               f' (R0 B Lamb ratio {B[0]/meas["Lamb"]:.1f}; R1 1S - 2S off {(R1[2] - meas["1S-2S"])/1e9:+.1f} GHz)')


def part_Ps():
    S = System(True)
    out.append('PART 3: positronium without spin, n = 2 (energies in alpha^2 Ha = alpha^4 m c^2 above -1/16)')
    ferrell = lambda n, L: 11/(64*n**4) - 1/(2*n**3*(2*L))
    d = {}; ok = True
    for L, lab in ((0.5, '2S'), (1.5, '2P')):
        Er = solve(S, L, lambda q: 2*q['D'] - 2*q['Wr'], -1/16)
        Eb = solve(S, L, lambda q: q['Jr'] - (2 - L), -1/16)
        a = lambda E: (E + 1/16)/ALPHA**2
        q = S.period(Eb, L); dev = (a(Eb) - ferrell(2, L))/abs(ferrell(2, L)); ok &= abs(dev) < 1e-3
        out.append(f'  {lab} (L {L}): E_res {a(Er):+.8f} | EBK {a(Eb):+.8f} vs Ferrell {ferrell(2, L):+.8f} (dev {dev:+.1e}) | at EBK: 2D/Omega_r {2*q["D"]/q["Wr"]:.9f},'
                   f' -E/D {-Eb/q["D"]:.9f}')
        d[lab] = (Er - Eb)/ALPHA**2
    shift = (d['2S'] - d['2P'])*A4
    out.append(f'  EBK = Ferrell within 1e-3: {ok}; delta(2S) - delta(2P) = {units(shift)} (prediction alpha^2/24 Ha = {units(A4/24)})')
    eng = {0: 18230.08e6, 1: 12744.86e6, 2: 8359.22e6, 'S1P1': 10917.86e6}
    meas = {0: (18501.02e6, 'Gurung et al. 2020'), 1: (13012.42e6, 'Hagena et al. 1993'), 2: (8627.94e6, 'Sheldon et al. 2023'), 'S1P1': (11180e6, 'Ley et al. 1994')}
    okD = True
    for k in (0, 1, 2, 'S1P1'):
        new = eng[k] + shift; m, src = meas[k]; okD &= abs(new - m) > 1e9
        lab = f'o-Ps 2^3S1 - 2^3P{k}' if k != 'S1P1' else 'p-Ps line 2^3S1 - 2^1P1'
        out.append(f'  {lab:24s}: engine {eng[k]/1e6:,.2f} -> with the joint-clock resonance {new/1e6:,.2f} MHz vs measured {m/1e6:,.2f} ({src}): {units(new - m)}')
    out.append(f'H717d (every Ps interval > 1 GHz from measurement): {okD}')


if __name__ == '__main__':
    part_H(); part_fs(); part_Ps()
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana717.txt', 'w').write(text + '\n')

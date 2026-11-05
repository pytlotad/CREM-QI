# Audit 713 (offline): the Darwin term from the geometry of the charge's fast motion (loop of the magnetic moment, zitterbewegung) against the
# measured Lamb shift.  Atomic units (m_e = hbar = e = a0 = 1, Ha), hydrogen with reduced mass mu in the density, positronium mu = 1/2.
# Geometries: G1 loop R = (g/2) hbar/mc, G2 loop r* = R/2, G3 isotropic shell (sqrt 3/2) hbar/mc, G4 the model's Plummer eps on the charge.
# Readings: (D) density |psi|^2 (quantum spin 1/2: rank-2 spin tensor = 0); (K) model: Langer orbit + classical spin with rule-Q lengths.
# usage: python3 ana713.py
import os, math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ALPHA = 7.2973525693e-3; MP = 1836.15267343; G2E = 1.00115965218128   # g/2 of the electron (= positron)
HA_HZ = 6.579683920502e15; HZ_EV = 4.135667696e-15; C_SI = 299792458.0
LAMB = 1057.845e6; LAMB_B = 1057.8298e6; FS_MEAS = 9911.200e6 + 1057.845e6   # Lundeen-Pipkin 1981; Bezginov 2019; Hagley-Pipkin 1994 + LP
LAMBDA = ALPHA                       # hbar/(m c) in a0
R1 = G2E*LAMBDA; RSTAR = R1/2; S3 = math.sqrt(3)/2*LAMBDA; EPS = 0.96682*RSTAR
out = []


def units(ha):
    f = ha*HA_HZ
    lam = (f'{C_SI/abs(f)*100:.4g} cm' if f != 0 else '-')
    return f'{f/1e6:+.2f} MHz = {f*HZ_EV*1e6:+.4f} ueV ({lam})'


def rho_S(n, mu):
    if n == 1: return lambda r: mu**3/math.pi*math.exp(-2*mu*r)
    if n == 2: return lambda r: mu**3/(8*math.pi)*(1 - mu*r/2)**2*math.exp(-mu*r)
    raise ValueError


def R21sq(mu): return lambda r: mu**3/24*(mu*r)**2*math.exp(-mu*r)


def shell_shift(rho, s):
    """Phi(0) - Phi(s) for a spherical density: exact energy shift of a charge spread on radius s (S states, any orientation)"""
    return quad(lambda r: rho(r)*4*math.pi*r*r*(1/r - 1/s), 0, s, epsabs=0, epsrel=1e-13, limit=200)[0]


def plummer_shift(rho):
    """<1/r - 1/sqrt(r^2 + eps^2)> for a spherical density (radial weight 4 pi r^2 rho)"""
    f = lambda r: rho(r)*4*math.pi*r*r*(1/r - 1/math.sqrt(r*r + EPS*EPS)) if r > 0 else rho(0)*4*math.pi*0.0
    return sum(quad(f, a, b, epsabs=0, epsrel=1e-12, limit=400)[0] for a, b in ((0, 10*EPS), (10*EPS, 1.0), (1.0, 80.0)))


def E_D(n, mu, two=False):   # Dirac Darwin (contact), hydrogen electron only; Ps both particles
    return (ALPHA**2*mu**3/n**3) if two else ALPHA**2*mu**3/(2*n**3)


def part_D():
    mu = 1/(1 + 1/MP); ok_a = True
    out.append('PART D (density reading), hydrogen, reduced mass mu = %.9f; Dirac Darwin E_D(1S) = %s, E_D(2S) = %s'
               % (mu, units(E_D(1, mu)), units(E_D(2, mu))))
    res = {}
    for tag, s, m2 in (('G1 loop R = (g/2) hbar/mc', R1, R1**2), ('G2 loop r* = R/2', RSTAR, RSTAR**2), ('G3 shell (sqrt3/2) hbar/mc', S3, S3**2)):
        for n in (1, 2):
            ex = shell_shift(rho_S(n, mu), s); contact = 2*math.pi/3*m2*rho_S(n, mu)(0)*(1 - mu*s)
            dev = ex/contact - 1; ok_a &= abs(dev) < 1e-4
            out.append(f'  {tag:28s} {n}S: exact shift {units(ex)}, /E_D {ex/E_D(n, mu):.6f}, contact x (1 - mu s) deviation {dev:+.2e}')
            res[(tag[:2], n)] = ex
        p = quad(lambda r: R21sq(mu)(r)*r*r*(1/r - 1/s), 0, s, epsabs=0, epsrel=1e-12)[0]   # 2P monopole (quantum spin: tensor = 0)
        res[(tag[:2], 'P')] = p
        out.append(f'  {tag:28s} 2P monopole shift {units(p)} (tensor part 0 for spin 1/2)')
    out.append(f'H713a (exact = (2 pi/3) <delta^2> |psi(0)|^2 (1 - mu s) within 1e-4, G1-G3, 1S and 2S): {ok_a}')
    lam = {}
    for g in ('G1', 'G2', 'G3'):
        lam[g] = res[(g, 2)] - E_D(2, mu) - res[(g, 'P')]
        out.append(f'  {g}: 2S1/2 - 2P1/2 = dE(2S) - E_D(2S) - dE(2P1/2) = {units(lam[g])}; /measured Lamb 1057.845 MHz: {lam[g]*HA_HZ/LAMB:+.3f}')
    okb = 6.5e9 <= lam['G1']*HA_HZ <= 7.5e9 and lam['G1']*HA_HZ > 5*LAMB
    okc = -15e9 <= lam['G2']*HA_HZ <= -14e9
    okd = lam['G3'] < 0 and abs(lam['G3']*HA_HZ) < 200e6
    # G1 with a classical spin and the QM triangle (|L| = sqrt2, |S| = sqrt3/2) on the density's <r^-3> = mu^3/24, planar angular average
    for j, cb in (('1/2', -2/math.sqrt(6)), ('3/2', 1/math.sqrt(6))):
        P2 = 1.5*cb*cb - 0.5; t = -(R1**2/4)*(mu**3/24)*P2
        out.append(f'  G1 hybrid (density <r^-3>, classical spin, QM triangle cos(L,S) {cb:+.4f}) 2P{j}: tensor shift {units(t)}')
    # G4: Plummer on the charge
    p4s = plummer_shift(rho_S(2, mu)); p4p = plummer_shift(lambda r: R21sq(mu)(r)/(4*math.pi))
    lam4 = p4s - E_D(2, mu) - p4p
    out.append(f'  G4 Plummer eps = {EPS/LAMBDA:.5f} hbar/mc on the charge: dE(2S) {units(p4s)}, dE(2P) {units(p4p)}; 2S1/2 - 2P1/2 {units(lam4)}'
               f' ({lam4*HA_HZ/LAMB:.1f}x measured)')
    okf = lam4*HA_HZ > 50e9
    # radius of a loop that gives exactly the measured Lamb shift in reading D
    f = lambda Rr: shell_shift(rho_S(2, mu), Rr) - quad(lambda r: R21sq(mu)(r)*r*r*(1/r - 1/Rr), 0, Rr, epsabs=0, epsrel=1e-12)[0] - E_D(2, mu) - LAMB/HA_HZ
    Rfit = brentq(f, 0.5*LAMBDA, 1.2*LAMBDA, xtol=1e-16)
    fB = lambda Rr: f(Rr) + (LAMB - LAMB_B)/HA_HZ
    RfitB = brentq(fB, 0.5*LAMBDA, 1.2*LAMBDA, xtol=1e-16)
    out.append(f'  loop radius reproducing the measured Lamb shift (D): R_fit = {Rfit/LAMBDA:.8f} hbar/mc (Lundeen-Pipkin), {RfitB/LAMBDA:.8f} (Bezginov);'
               f' R_fit^2 = {(Rfit/LAMBDA)**2:.5f} vs Dirac 3/4 (excess {(Rfit/LAMBDA)**2/0.75 - 1:+.4f}); 1S at R_fit: dE - E_D {units(shell_shift(rho_S(1, mu), Rfit) - E_D(1, mu))}')
    out.append(f'H713b (G1: 2S1/2 - 2P1/2 in [+6.5, +7.5] GHz and > 5x measured): {okb}')
    out.append(f'H713c (G2: 2S1/2 - 2P1/2 in [-15, -14] GHz): {okc}')
    out.append(f'H713d (G3: 2S1/2 - 2P1/2 < 0, |.| < 200 MHz): {okd}')
    out.append(f'H713f (G4: 2S1/2 - 2P1/2 > +50 GHz): {okf}')
    # proton size, for scale
    rp = 0.8409e-15/5.29177210903e-11
    out.append(f'  scale: proton charge radius 0.8409 fm on 2S: {units(2*math.pi/3*rp**2*rho_S(2, mu)(0))}')


def kepler_ring_average(n, L, mu, R, cosb, M=6000, K=96):
    """time average over a Kepler orbit (reduced mass mu, k = 1) of the ring-averaged Coulomb potential minus the point one; ring of radius R in the
    plane perpendicular to the spin s, s at angle beta from L = z (azimuth 0.3 rad); exact numerical average over K ring points"""
    a = n*n/mu; e = math.sqrt(1 - (L/n)**2)
    u = (np.arange(M) + 0.5)*2*math.pi/M; r = a*(1 - e*np.cos(u))
    x = a*(np.cos(u) - e); y = a*math.sqrt(1 - e*e)*np.sin(u); w = (1 - e*np.cos(u))/M
    sb = math.sqrt(1 - cosb*cosb); s = np.array([sb*math.cos(0.3), sb*math.sin(0.3), cosb])
    e1 = np.cross(s, [1.0, 0, 0]); e1 /= np.linalg.norm(e1); e2 = np.cross(s, e1)
    ph = (np.arange(K) + 0.5)*2*math.pi/K
    d = R*(np.outer(np.cos(ph), e1) + np.outer(np.sin(ph), e2))            # K x 3
    P = np.stack([x, y, 0*x], 1)                                             # M x 3
    dist = np.linalg.norm(P[:, None, :] + d[None, :, :], axis=2)
    dV = -(1/dist).mean(1) + 1/r
    rp = a*(1 - e)
    return float(np.sum(dV*w)), rp


def part_K():
    mu = 1/(1 + 1/MP); ok = True; sh = {}
    out.append('PART K (model reading): Langer orbits, classical spin with rule-Q lengths (|s| = hbar, k = |L + s| = j + 1/2), loop G1 R = (g/2) hbar/mc')
    for name, n, L, k in (('1S band', 1, 0.5, 1), ('2S band', 2, 0.5, 1), ('2P1/2', 2, 1.5, 1), ('2P3/2', 2, 1.5, 2), ('3P1/2', 3, 1.5, 1),
                          ('3P3/2', 3, 1.5, 2), ('4S band', 4, 0.5, 1), ('4P1/2', 4, 1.5, 1)):
        cb = (k*k - L*L - 1)/(2*L); P2 = 1.5*cb*cb - 0.5
        num, rp = kepler_ring_average(n, L, mu, R1, cb)
        formula = -(R1**2/4)*mu**3/(n**3*L**3)*P2
        ok &= rp/R1 > 10
        sh[name] = num
        out.append(f'  {name:8s} n {n} L {L}: cos(s, L) {cb:+.4f}, periapsis/R {rp/R1:7.1f} (contact 0); ring average {units(num)}'
                   f' | quadrupole formula -(R^2/4) <r^-3> P2 {units(formula)} (ratio {num/formula:.5f})')
    dfs = sh['2P3/2'] - sh['2P1/2']; ds = sh['2S band'] - sh['2P1/2']
    out.append(f'  2P3/2 - 2P1/2 change {units(dfs)} = {dfs*HA_HZ/FS_MEAS*100:+.1f} % of measured 10969.045 MHz;'
               f' 2S1/2 - 2P1/2 change (on top of import 21) {units(ds)} = {ds*HA_HZ/LAMB:+.1f}x measured Lamb')
    okE = ok and dfs*HA_HZ > 2e9 and sh['2S band']*HA_HZ > 30e9
    out.append(f'H713e (contact 0 on all orbits, 2P3/2 - 2P1/2 change > +2 GHz, 2S band > +30 GHz): {okE}')


def part_Ps():
    mu = 0.5
    out.append('PART Ps (density reading, G1 loops on both particles, coplanar rings -- o-Ps spins parallel, p-Ps antiparallel -- independent phases:'
               ' relative smear d = 2R|sin(u/2)|)')
    for n in (1, 2):
        rho = rho_S(n, mu)
        ex = quad(lambda u: shell_shift(rho, 2*R1*abs(math.sin(u/2))) if abs(math.sin(u/2)) > 1e-12 else 0.0, 0, 2*math.pi, epsabs=0, epsrel=1e-11,
                  limit=200)[0]/(2*math.pi)
        ed = E_D(n, mu, two=True)
        for ch, st in (('p-Ps', f'{n}^1S0'), ('o-Ps', f'{n}^3S1')):
            out.append(f'  {ch} {st}: smear shift {units(ex)}; Darwin of both particles {units(ed)}; excess {units(ex - ed)} (/E_D {ex/ed:.5f})')
    out.append('  (contact independent of spin: the two channels equal by construction; the 2P part in reading K needs the channel spin configuration)')


if __name__ == '__main__':
    part_D(); part_K(); part_Ps()
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana713.txt', 'w').write(text + '\n')

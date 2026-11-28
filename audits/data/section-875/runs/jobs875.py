# Audit 875: job list (one line per run: name then env assignments), TEND in orbits T = 2 pi N0^3 / mu.
import math
MH = 1836.15267343; PS = 'M1=1 G1=2.00231930436'
CS = {'10': 10.0, '30': 30.0, '137': 137.035999084}
def T(m1, n0): mu = m1 / (m1 + 1); return 2 * math.pi * n0 ** 3 / mu
L = []
def job(name, sp, n, c, norb, extra=''):
    m1 = 1 if sp == 'P' else MH; n0 = n + 0.02
    L.append(f"{name} MODE=1 WDELAY=1 LL=1 RCUT=0.1 C={CS[c]} {PS + ' ' if sp == 'P' else ''}N0={n0:.2f} E0=0 TEND={norb * T(m1, n0):.0f} {extra}".strip())
for sp in ('H', 'P'):
    for n in (1, 2):
        for c in ('10', '30', '137'):
            job(f"F_{sp}{n}_c{c}", sp, n, c, 6000); job(f"K_{sp}{n}_c{c}", sp, n, c, 1000, 'KAPPA=0')
for n in (1, 2):
    job(f"S_H{n}_F0", 'H', n, '137', 6000, 'MS1=-0.5 MS2=0.5'); job(f"S_H{n}_F1", 'H', n, '137', 6000, 'MS1=0.5 MS2=0.5')
    job(f"S_P{n}_pPs", 'P', n, '137', 6000, 'MS1=-0.5 MS2=0.5'); job(f"S_P{n}_oPs", 'P', n, '137', 6000, 'MS1=0.5 MS2=0.5')
L.sort(key=lambda l: -float(l.split('TEND=')[1].split()[0]))
print('\n'.join(L))

"""Audit 628, uzupelnienie 2: split of the EBK alpha^4 shift into the kinetic part (Darwin off, DARWIN=0) and the orbit-orbit part
(EBK with Darwin - EBK without), against QM kinetic (two particles, p^4) and Ferrell - QM kinetic.  usage: python3 decomp628.py out_dir"""
import sys


def load(path):
    c, L = {}, {}
    for line in open(path):
        w = line.split()
        if w[0] == 'CONST':
            c = {w[i]: float(w[i + 1]) for i in range(1, len(w), 2)}
        if w[0] == 'LEVEL':
            L[(int(w[2]), int(w[4]))] = float(w[10])
    return c, L


d = sys.argv[1]
for s, full in (('0.5', 'levels_ps_s0.5.txt'), ('1', 'levels_ps.txt')):
    c, nod = load('%s/lv_nod_s%s.txt' % (d, s))
    _, wd = load('%s/%s' % (d, full))
    m, k, hb, cc, mu = c['m1'], c['k'], c['hbar'], c['c'], c['mu']
    a, al = hb * hb / (mu * k), k / (hb * cc)
    for (n, l), dE in sorted(nod.items()):
        E0 = -mu * k * k / (2 * hb * hb * n * n)
        r1, r2 = 1 / (n * n * a), 1 / (n ** 3 * (l + 0.5) * a * a)
        kin = -(2 / m ** 3) * 4 * mu * mu * (E0 * E0 + 2 * E0 * k * r1 + k * k * r2) / (8 * cc * cc)
        fer = m * al ** 4 * cc * cc * (11 / (64 * n ** 4) - 1 / (2 * n ** 3 * (2 * l + 1)))
        print("s=%s n=%d l=%d  EBK_bez_Darwina/kin_QM - 1 = %+.3e   (EBK_z - EBK_bez)/(Ferrell - kin_QM) - 1 = %+.3e"
              % (s, n, l, dE / kin - 1, (wd[(n, l)] - dE) / (fer - kin) - 1))

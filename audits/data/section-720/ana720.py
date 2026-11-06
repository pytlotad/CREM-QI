"""Audit 720: the spin-independent m alpha^5 block of Ps (Pineda & Soto, PRD 59, 016005 (1998), eq. (3.8); Fulton & Martin 1954) as an import,
with Bethe logarithms ln k0(n, l) from Jentschura & Mohr, PRA 72, 012110 (2005) (arXiv quant-ph/0504002; = Drake & Swainson, PRA 41, 1243 (1990)
for n <= 20), and the S-state convention shift (3/4) T_n that brings the engine's one-loop HFS split (575: triplet +T/4, singlet -3T/4) to the formula's
(triplet T, singlet 0).  Import per level:
    Delta(n, 0) = P { -(2 alpha/3 pi)/n^3 [9 ln alpha + 7 ln n + 8 ln k0 - 14 ln 2 - 49/15 - 7 (H_n + (n - 1)/(2n))] + (3/2)(alpha/pi)(-32/9 - 2 ln 2)/n^3 }
    Delta(n, l >= 1) = P { -(16 alpha/3 pi)/n^3 [ln k0 + 7/(16 l (l + 1)(2l + 1))] },   P = m c^2 alpha^4/8.
usage: python3 ana720.py [p640 output of the new engine] -> predictions (always) and engine vs formula (with the file)."""
import math, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
a = 7.2973525693e-3; mc2 = 1.2355899662e20/1e6; P = mc2*a**4/8; ap = a/math.pi   # MHz
lnk0 = {}
for line in open(os.path.join(HERE, 'lnk0_n30.tsv')):
    n, l, v = line.split(); lnk0[(int(n), int(l))] = float(v)


def L(n, l):
    nn = min(n, 30); return lnk0[(nn, min(l, nn - 1))]


def H(n): return sum(1/k for k in range(1, n + 1))


def SI(n, l):
    if l == 0:
        return P*(-(2*a/(3*math.pi))/n**3*(9*math.log(a) + 7*math.log(n) + 8*L(n, 0) - 14*math.log(2) - 49/15 - 7*(H(n) + (n - 1)/(2*n))))
    return P*(-(16*a/(3*math.pi))/n**3*(L(n, l) + 7/(16*l*(l + 1)*(2*l + 1))))


def T(n): return P*(14/3)/n**3*(3*a/(7*math.pi))*(-32/9 - 2*math.log(2))


def imp(n, l): return SI(n, l) + (0.75*T(n) if l == 0 else 0.0)


MHz_per_ueV = 241.79892623   # 1 ueV = 241.798926 MHz
out = []
out.append('P = m c^2 alpha^4/8 = %.4f MHz' % P)
for n, l in [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2), (4, 0), (10, 0), (30, 0), (100, 0)]:
    out.append('import n %3d l %d: SI %+12.4f  conv %+10.4f  total %+12.4f MHz (%+.5f ueV)' % (n, l, SI(n, l), imp(n, l) - SI(n, l), imp(n, l), imp(n, l)/MHz_per_ueV))
dS2P = imp(2, 0) - imp(2, 1); d12 = imp(2, 0) - imp(1, 0)
out.append('2S - 2P import (adds to every nu_J and to 2^3S1 - 2^1P1): %+.4f MHz' % dS2P)
out.append('2S - 1S import (adds to 1S - 2S): %+.4f MHz (629 missing alpha^5 without HFS part: -1721.4)' % d12)
eng = {'nu0': 18230.08, 'nu1': 12744.86, 'nu2': 8359.22, 'SP': 10917.86}
meas = {'nu0': (18501.02, 0.61, 'Gurung et al., PRL 125, 073002 (2020)'), 'nu1': (13012.42, 1.67, 'Hagena et al., PRL 71, 2887 (1993)'),
        'nu2': (8627.94, 0.95, 'Sheldon et al., PRL 131, 043001 (2023)'), 'SP': (11180.0, None, 'Ley et al., Hyperfine Interact. 89, 327 (1994); uncertainty not available here (only Fig. 5 of Karshenboim 2004)')}
for k in ('nu0', 'nu1', 'nu2', 'SP'):
    new = eng[k] + dS2P; m, u, ref = meas[k]
    us = '(%.2f)' % u if u else ''; sg = '; %.1f sigma' % (abs(new - m)/u) if u else ''
    out.append('%s: engine (714) %.2f -> predicted %.2f MHz vs %.2f%s [%s]: residual %+.2f -> %+.2f MHz (%+.4f ueV%s)'
               % (k, eng[k], new, m, us, ref, eng[k] - m, new - m, (new - m)/MHz_per_ueV, sg))
r12 = 1737.0 + d12
out.append('1S - 2S: residual (629/638) +1737.0 -> predicted %+.2f MHz vs Fee et al. 1993 (1 233 607 216.4(3.2) MHz)' % r12)
# 2p -> 1s line shift (o-Ps and p-Ps): import(2P) - import(1S)
dl = imp(2, 1) - imp(1, 0)
out.append('2p -> 1s line shift: %+.4f MHz = %+.5f ueV' % (dl, dl/MHz_per_ueV))
if len(sys.argv) > 1:
    new = {}
    for line in open(sys.argv[1]):
        w = line.split(); new[w[0]] = float(w[1])
    old = {}
    for line in open(os.path.join(HERE, '../section-714/out/p640_714.txt')):
        w = line.split(); old[w[0]] = float(w[1])
    ok = True
    for k, (n, l) in {'2^3S1': (2, 0), '2^1P1': (2, 1), '2^3P0': (2, 1), '2^3P1': (2, 1), '2^3P2': (2, 1)}.items():
        d = new[k] - old[k]; ok &= abs(d - imp(n, l)) < 0.01
        out.append('engine %s: new - old %+.4f MHz vs import %+.4f (diff %+.1e)' % (k, d, imp(n, l), d - imp(n, l)))
    out.append('H720a (engine level change = import within 0.01 MHz, n = 2): %s' % ok)
text = '\n'.join(out); print(text)
os.makedirs(os.path.join(HERE, 'out'), exist_ok=True)
open(os.path.join(HERE, 'out', 'ana720.txt' if len(sys.argv) > 1 else 'pred720.txt'), 'w').write(text + '\n')

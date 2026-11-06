"""Audit 724: engine constants CODATA 2018 -> 2022 (epsilon0, hbar = h/2pi, bohrRadius).  Reads the probes before/after (out/p721_*, p640_*, p629_*)
and prints the engine's own intervals and residuals against measurements.  usage: python3 ana724.py -> out/ana724.txt"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
def rd(fn):
    d = {}
    for line in open(os.path.join(HERE, 'out', fn)):
        w = line.split()
        if len(w) >= 2:
            try: d[w[0]] = float(w[1])
            except ValueError: pass
    return d
uev = 241.79892623; c = 299792458.0; out = []
Rc22 = 10973731.568157*c/1e6
for tag in ('before', 'after'):
    E = rd('p721_%s.txt' % tag); a = 1/E['alpha_inv']; mc2 = E['me_c2_MHz']; x = E['mp_over_me']/(1 + E['mp_over_me'])
    R = a*a*mc2/2
    lamb = E['F_2S'] - E['F_2P1/2']; s3 = E['F_2P3/2'] - E['F_2S']; fs = E['F_2P3/2'] - E['F_2P1/2']
    c12 = 0.75*E['binding_MHz'] + E['F_2S'] - E['F_1S']
    out.append('%s: alpha^-1 %.9f  R_inf c %.4f MHz (CODATA 2022 %.4f: %+.4f MHz, %+.2e)  binding %.6f MHz' % (tag, 1/a, R, Rc22, R - Rc22, (R - Rc22)/Rc22, E['binding_MHz']))
    out.append('  H 2S1/2 - 2P1/2 %.4f MHz (Bezginov 2019 1057.8298(32): %+.4f);  2P3/2 - 2S1/2 %.4f (Hagley-Pipkin 9911.200(12): %+.4f);  2P3/2 - 2P1/2 %.4f (10969.045(15): %+.4f)'
               % (lamb, lamb - 1057.8298, s3, s3 - 9911.200, fs, fs - 10969.045))
    out.append('  H 1S - 2S centroid %.4f MHz vs 2466061413.187035(10) (Parthey 2011): %+.4f MHz (%+.4f ueV)' % (c12, c12 - 2466061413.187035, (c12 - 2466061413.187035)/uev))
    # line 2p -> 1s (spin centroid of each j, F = 1 hyperfine as in [W]: level differences without U on 2P)
    for j, k in (('1/2', 'F_2P1/2'), ('3/2', 'F_2P3/2')):
        f = 0.75*E['binding_MHz'] + E[k] - E['F_1S']
        out.append('  H 2P%s - 1S (no hyperfine) %.4f MHz = %.9f eV = %.6f nm' % (j, f, f/uev*1e-6, c/(f*1e6)*1e9))
    P = rd('p640_%s.txt' % tag)
    nu = [P['2^3S1'] - P['2^3P%d' % J] for J in range(3)]; sp = P['2^3S1'] - P['2^1P1']
    out.append('  Ps (o-Ps) nu_0,1,2 %.4f / %.4f / %.4f MHz (Gurung 18501.02(61): %+.4f; Hagena 13012.42(1.67): %+.4f; Sheldon 8627.94(95): %+.4f)'
               % (nu[0], nu[1], nu[2], nu[0] - 18501.02, nu[1] - 13012.42, nu[2] - 8627.94))
    out.append('  Ps (p-Ps) 2^3S1 - 2^1P1 %.4f MHz (Ley 1994 11180: %+.4f)' % (sp, sp - 11180.0))
    t = open(os.path.join(HERE, 'out', 'p629_%s.txt' % tag)).read().split('total')[1].split()[0]
    ps12 = float(t)
    out.append('  Ps (o-Ps) 1^3S1 - 2^3S1 %.4f MHz vs 1233607216.4(3.2) (Fee 1993): %+.4f MHz (%+.4f ueV)' % (ps12, ps12 - 1233607216.4, (ps12 - 1233607216.4)/uev))
text = '\n'.join(out); print(text)
open(os.path.join(HERE, 'out', 'ana724.txt'), 'w').write(text + '\n')

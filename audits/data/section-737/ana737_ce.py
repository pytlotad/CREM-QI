"""Audit 737, parts C' and E (run separately after the main run stopped on an orbital-label IndexError at n = 10 in part C'; parts A, D, B of
ana737.py are complete in out/ana737.log).  Same functions as ana737.py."""
import os, time
import ana737 as m
P = m.P; t0 = time.time()
TR = [((2, 1), (1, 0)), ((3, 1), (1, 0)), ((3, 1), (2, 0))]
nm = lambda i, f: '%d%s -> %d%s' % (i[0], 'spdfghiklmnoqrtuvwxyz'[i[1]], f[0], 'spdfghiklmnoqrtuvwxyz'[f[1]])
P('\n(C\') POST HOC replacement of H737c: order-1 gain vs n -- circular n -> n-1 and np -> (n-1)s, |R/R_QM - 1| at k0 and k1-in')
for n in (5, 10, 20):
    for (i, f) in (((n, n - 1), (n - 1, n - 2)), ((n, 1), (n - 1, 0))):
        q, e0, _, _ = m.radial(*i, *f, 0); _, e1, _, _ = m.radial(*i, *f, 1)
        P('  %s: k0 %.2e  k1 %.2e' % (nm(i, f), abs(e0/q - 1), abs(e1/q - 1)))
    m._cache.clear()
P('\n(E) grid 2x finer, k = 1 (in, gl); base values recomputed here')
ok = True
for (i, f) in TR:
    m._cache.clear()
    q, e1, _, _ = m.radial(*i, *f, 1); _, eg, _, _ = m.radial(*i, *f, 1, glob=True); b = (e1/q, eg/q)
    m._cache.clear()
    q, e1, _, _ = m.radial(*i, *f, 1, fine=2); _, eg, _, _ = m.radial(*i, *f, 1, fine=2, glob=True)
    d = (abs(e1/q - b[0]), abs(eg/q - b[1])); ok &= max(d) < 1e-5
    P('  %s: base in %.7f gl %.7f; fine in %.7f gl %.7f; diff %.1e %.1e' % (nm(i, f), *b, e1/q, eg/q, *d))
P('H737e (order 0 = 733 to 1e-5 [part A: True] and grid 2x diff at k=1 < 1e-5): %s' % ok)
P('\nruntime %.0f s' % (time.time() - t0))
open(os.path.join(m.HERE, 'out', 'ana737_ce.txt'), 'w').write('\n'.join(m.out) + '\n')

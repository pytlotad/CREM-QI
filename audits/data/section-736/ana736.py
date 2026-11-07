"""Audit 736, POST HOC: the E2 channel 3d -> 2s left out of H736d (QM vs EBK, same functions as ana735) and the proton term of the 21 cm M1
check (engine sum (gamma_e - gamma_p)^2 hbar^2/4 vs electron only in 735)."""
import os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, '../section-735/ana735.py')).read().split('for nC in (10, 20, 40)')[0]   # functions only (its own header lines print once)
g = {'__file__': os.path.join(HERE, '../section-735/ana735.py')}; exec(src, g)
out = []
def P(s): out.append(s); print(s, flush=True)
for (ni, li, nf, lf) in ((3, 2, 1, 0), (3, 2, 2, 0)):
    q, e = g['rad'](ni, li, nf, lf, 2); G = g['ang_E2'](li, lf)
    P('%d%d -> %d%d: R2 QM %+.5f EBK %+.5f ratio %.5f; A_E2 QM %.4e EBK %.4e s^-1 (G %.5f)' % (ni, li, nf, lf, q, e, e/q, g['A_E2'](ni, li, nf, lf, q), g['A_E2'](ni, li, nf, lf, e), G))
q1, _ = g['rad'](3, 2, 1, 0, 2); q2, _ = g['rad'](3, 2, 2, 0, 2); r1, _ = g['rad'](3, 2, 2, 1, 1)
tq = g['A_E2'](3, 2, 1, 0, q1) + g['A_E2'](3, 2, 2, 0, q2); e1 = g['A_E1'](3, 2, 2, 1, r1)
P('3d total E2 (QM) %.4e s^-1, A_E1(3d -> 2p) QM %.4e -> E2 share of the count %.4e (735 quoted 3d -> 1s only: 9.164e-6)' % (tq, e1, tq/e1))
ge, gp = 1.76085963023e11, 2.6752218744e8
P('21 cm: (gamma_e - gamma_p)^2/gamma_e^2 - 1 = %.4e (proton moment in the M1 element; 735 reference electron only)' % ((1 + gp/ge)**2 - 1))
open(os.path.join(HERE, 'out', 'ana736.txt'), 'w').write('\n'.join(out) + '\n')

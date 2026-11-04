# Audit 704 POST HOC: the o-Ps tensor energy alone (the p-Ps classical U_dd of parallel moments, absent from the quantum singlet, not subtracted)
# against the measured b S12 at the classical vectors; expectation 8/9 (classical <r^-3> at L = 1.5 hbar).
# usage: python3 post704.py <702 runs dir with s8192/>
import sys, os
import numpy as np
from ana704 import parse, B_MHZ, HEV, EV
R2 = sys.argv[1]; out = []
for s in ('0', '1', '2', '3'):
    o = parse(f'{R2}/s8192/ps_A2_{s}.err'); p = parse(f'{R2}/s8192/ps_A1_{s}.err')
    L, s1, s2 = o['L'], o['s1'], o['s2']; S12 = 0.5*float(s1 @ s2) - 1.5*float(s1 @ L)*float(s2 @ L)
    Do = np.mean([q['Udd'] for q in o['k']]); Dp = np.mean([q['Udd'] for q in p['k']]); Eb = B_MHZ*1e6*HEV*S12*EV
    out.append(f'Ps seed {s}: o-Ps <U_dd> {Do/EV*1e6:+.4f} ueV / measured b S12 {Eb/EV*1e6:+.4f} ueV = {Do/Eb:.4f} (8/9 = {8/9:.4f}) | p-Ps classical <U_dd> {Dp/EV*1e6:+.4f} ueV (quantum singlet: 0)')
text = '\n'.join(out); print(text); os.makedirs('out', exist_ok=True); open('out/post704.txt', 'w').write(text + '\n')

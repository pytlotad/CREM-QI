# Audit 779, POST HOC (3b): phase error of the 2000-step scheme (same step as guide772/predict776) on the EBK tori, own mode, CAUST -pi/2 -> out/post779_ecc.txt
import os, sys, io, contextlib, numpy as np; D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
with contextlib.redirect_stdout(io.StringIO()): import ana779 as A9
from phase779 import wrap
for b in ('0.05', '0.2', '0.8'):
    R = [(k, h) for k, h in A9.E.items() if k.split(':')[1] == b and k.split(':')[2] != '0']
    bad = []
    for k, h in R:
        own = h['modes'].index((int(round(h['J'] - 0.5)), int(round(h['L'])))); e = abs(float(wrap(h['Pc'][own] - h['Pn'][own])))
        Lam = np.sqrt(h['L'] ** 2 + 2 * float(b)); ecc = np.sqrt(1 - (Lam / h['N']) ** 2)
        if e > 1e-3: bad.append((int(round(h['J'] - 0.5)), int(round(h['L'])), round(e, 4), round(ecc, 4)))
    print(f"beta {b}: torusy EBK z bledem fazy wlasnego modu > 1e-3: {len(bad)}/{len(R)}: {bad}")

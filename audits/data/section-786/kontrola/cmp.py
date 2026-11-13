# Consistency check secular vs mechanical on HEAD vs 699/700 (same quantities as ana700: secB0/mechB0 rate ratios L^, lepton, heavy; c = population drift ratio).
import sys, os, glob, numpy as np
sys.path.insert(0, '/home/teddy/Projekty/CREM-QI/audits/data/section-699')
from ana699 import parse, pops
S = os.path.dirname(os.path.abspath(__file__))
def q(f):
    d = parse(f); g = d['inc']['sec']; X0 = [g['L0'], g['s10'], g['s20']]
    r = [np.linalg.norm(np.asarray(d['inc']['secB0'][a]) - g[b]) / np.linalg.norm(np.asarray(d['inc']['mechB0'][a]) - g[b]) for a, b in (('L1', 'L0'), ('s21', 's20'), ('s11', 's10'))]
    w, V = np.linalg.eigh(d['H']['all']); p0 = pops(V, X0)
    dD = {k: pops(V, [d['inc'][k]['L1'], d['inc'][k]['s11'], d['inc'][k]['s21']]) - p0 for k in ('sec', 'mech')}
    c = float(dD['mech'] @ dD['sec'] / (dD['sec'] @ dD['sec']))
    return r, c, float(g['L0'] @ g['s20'])
for v, old in (('both', 'old700/both'), ('def', 'old699/def')):
    print(f"== {v} (HEAD wobec {old}) ==")
    for f in sorted(glob.glob(f"{S}/{v}/*.err")):
        b = os.path.basename(f); o = f"{S}/{old}/{b}"
        rn, cn, ls = q(f); ro, co, lo = q(o) if os.path.exists(o) else ([np.nan] * 3, np.nan, np.nan)
        print(f"  {b[:-4]:26s} L.s {ls:+.3f} ({lo:+.3f}) | sec/mech L^ {rn[0]:.4f} ({ro[0]:.4f}) lepton {rn[1]:.4f} ({ro[1]:.4f}) heavy {rn[2]:.4f} ({ro[2]:.4f}) | c {cn:.3f} ({co:.3f})")

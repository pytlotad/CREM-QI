"""Audit 632: where the few-MHz difference of the 2P triplet level between the basin centre (Born turn, 631) and positions inside the
basin comes from -- emission-ledger fields of the 3d -> 2p and 2p -> 1s photons, grouped by J' and by whether the spins were turned.
usage: python3 ledger632.py dir_with_new_*.err"""
import glob
import re
import sys

hbar = 1.054571817e-34
h = 2 * 3.141592653589793 * hbar
VEC = {'dir', 'sphoton', 'ax', 'peri', 'm1', 'm2', 'vcm', 'hidJ'}


def parse(line):
    w = line.split()
    rec, side, j = {}, '', 1
    while j < len(w):
        k = w[j]
        if k == '|':
            side = w[j + 1] + '_'
            j += 2
            continue
        n = 3 if k in VEC else 1
        rec[side + k] = [float(x) for x in w[j + 1:j + 1 + n]]
        j += 1 + n
    return rec


runs = []
for p in sorted(glob.glob(sys.argv[1] + '/new_*.err')):
    cur = None
    for line in open(p):
        if line.startswith('START idx'):
            cur = {'idx': int(line.split()[2]), 'ph': None, 'basin': None, 'J': None, 'led': []}
            runs.append(cur)
        elif cur and line.startswith('BORN formation'):
            cur['J'] = int(re.search(r"J=(\d+)", line).group(1))
        elif cur and line.startswith('BORN photon'):
            m = re.search(r"J=(\d+) -> l'=\d+ J'=(\d+) \(basin (\d+)", line)
            cur['ph'], cur['basin'] = int(m.group(2)), int(m.group(3))
        elif cur and line.startswith('CREM_LEDGER t'):
            cur['led'].append(parse(line))


def rng(sel, key, i, scale):
    v = [r['led'][i][key][0] * scale for r in sel]
    return "%.6f..%.6f" % (min(v), max(v))


for Jp in (0, 1, 2):
    for turned in (False, True):
        sel = [r for r in runs if r['ph'] == Jp and (r['ph'] != r['basin']) == turned and len(r['led']) > 1]
        if not sel:
            continue
        print("J'=%d obrot=%s n=%d (J przy powstaniu %s)" % (Jp, turned, len(sel), sorted(set(r['J'] for r in sel))))
        print("   po 3d->2p:    L/hbar %s   UF %s MHz" % (rng(sel, 'A_L', 0, 1 / hbar), rng(sel, 'A_UF', 0, 1e-6 / h)))
        print("   przed 2p->1s: L/hbar %s   UF %s MHz   E %s MHz" % (rng(sel, 'B_L', 1, 1 / hbar), rng(sel, 'B_UF', 1, 1e-6 / h),
                                                              rng(sel, 'B_E', 1, 1e-6 / h)))
        print("   po 2p->1s:    UF %s MHz   E %s MHz   Eg %s MHz" % (rng(sel, 'A_UF', 1, 1e-6 / h), rng(sel, 'A_E', 1, 1e-6 / h),
                                                            rng(sel, 'Eg', 1, 1e-6 / h)))

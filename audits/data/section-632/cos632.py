"""Audit 632, H632b-c: cos(s1, s2) of the two spins (from the moments in the emission ledger, s_i = sign(q_i) m_i) before / after each photon,
and the 1S level energy U + F after the 2p -> 1s photon against it.  usage: python3 cos632.py dir_with_new_*.err"""
import glob
import math
import re
import sys

h = 2 * math.pi * 1.054571817e-34
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


def cos_spins(m1, m2):
    # electron (first, q < 0): s1 = -m1; positron (second): s2 = +m2
    s1 = [-x for x in m1]
    dot = sum(a * b for a, b in zip(s1, m2))
    return dot / (math.sqrt(sum(a * a for a in s1)) * math.sqrt(sum(b * b for b in m2)))


rows = []
for p in sorted(glob.glob(sys.argv[1] + '/new_*.err')):
    cur = None
    for line in open(p):
        if line.startswith('START idx'):
            cur = {'ph': None, 'basin': None, 'led': []}
            rows.append(cur)
        elif cur and line.startswith('BORN photon'):
            m = re.search(r"J'=(\d+) \(basin (\d+)", line)
            cur['ph'], cur['basin'] = int(m.group(1)), int(m.group(2))
        elif cur and line.startswith('CREM_LEDGER t'):
            cur['led'].append(parse(line))
rows = [r for r in rows if len(r['led']) > 1]
data = []
for r in rows:
    c = [cos_spins(r['led'][0]['B_m1'], r['led'][0]['B_m2']), cos_spins(r['led'][0]['A_m1'], r['led'][0]['A_m2']),
         cos_spins(r['led'][1]['B_m1'], r['led'][1]['B_m2']), cos_spins(r['led'][1]['A_m1'], r['led'][1]['A_m2'])]
    data.append((r, c, r['led'][1]['A_UF'][0] / h * 1e-6))
print("H632c: 1 - cos(s1, s2) [przed 3d->2p | po 3d->2p | przed 2p->1s | po 2p->1s]:")
for Jp in (0, 1, 2):
    for turned in (False, True):
        sel = [(r, c, u) for r, c, u in data if r['ph'] == Jp and (r['ph'] != r['basin']) == turned]
        if not sel:
            continue
        cols = [[1 - c[i] for _, c, _ in sel] for i in range(4)]
        print("  J'=%d obrot=%-5s n=%2d " % (Jp, turned, len(sel)) + " | ".join("%.2e..%.2e" % (min(x), max(x)) for x in cols)
              + "   UF_1S %.3f..%.3f MHz" % (min(u for *_, u in sel), max(u for *_, u in sel)))
# H632b: UF_1S against cos after the last photon, all runs
xs = [c[3] for _, c, _ in data]
ys = [u for *_, u in data]
n = len(xs)
mx, my = sum(xs) / n, sum(ys) / n
sxx = sum((x - mx) ** 2 for x in xs)
syy = sum((y - my) ** 2 for y in ys)
sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
slope = sxy / sxx if sxx > 0 else float('nan')
r = sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float('nan')
print("H632b: UF_1S wobec cos(s1, s2) po 2p->1s: n = %d, korelacja r = %.6f, nachylenie %.3f MHz na jednostke cos, wyraz wolny przy cos = 1: %.3f MHz"
      % (n, r, slope, my + slope * (1 - mx)))

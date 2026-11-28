# Audit 875: local field + Landau-Lifshitz + LW -- verdicts H875c, r, b, m, k, t, z, e -> out/ana875.txt.
# Written after prereg875.txt and the launch, before any result.  Extra diagnostic (not a rule, noted before results): for Ps the
# self term -W_rr is only ~half of the dipole power (the cross term is in the LW forces), so H875z as written compares the self part;
# the full balance W_g - dE per cycle is printed alongside.
import os, glob, math, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
RINF = 10973731.568157 * 299792458.0; C0 = 299792458.0; HEV = 4.135667696e-15; EH = 2 * RINF
def load(p):
    hdr = {}; rows = []; stop = fall = False
    for ln in open(p):
        if ln.startswith('# 871'):
            t = ln.split()[2:]; hdr = {t[i]: float(t[i + 1]) for i in range(0, len(t) - 1, 2)}
        elif ln.startswith('# STOP fall'): stop = fall = True
        elif ln.startswith('# STOP'): stop = True
        elif not ln.startswith('#') and ln.strip(): rows.append([float(x) for x in ln.split()])
    return hdr, (np.array(rows) if rows else np.zeros((0, 15))), stop, fall
def summ(p):
    name = os.path.basename(p)[:-4]; hdr, q, stop, fall = load(p); m = len(q)
    mu = hdr['M1'] * hdr['M2'] / (hdr['M1'] + hdr['M2']); n = round(hdr['N0']); c = hdr['C']; T = 2 * math.pi * hdr['N0'] ** 3 / mu
    s0 = -4 * math.pi / (3 * c ** 3 * n ** 2)
    ch = 'bez spinu' if hdr['MS1'] == 0 and hdr['MS2'] == 0 else {('H', -1): 'F=0', ('H', 1): 'F=1', ('Ps', -1): 'p-Ps', ('Ps', 1): 'o-Ps'}[('H' if hdr['M1'] > 2 else 'Ps', int(round(2 * hdr['MS1'])))]
    r = dict(name=name, sp='H' if hdr['M1'] > 2 else 'Ps', n=n, c=c, mu=mu, T=T, s0=s0, kappa=hdr['KAPPA'], cyc=m, stop=stop, fall=fall, ch=ch,
             orb=(q[-1, 1] / T if m else 0))
    if m < 10: r.update(ok=False, cls='PEKA' if stop else '?', s=float('nan'), J=float('nan')); return r
    sel = slice(0, m) if stop else slice(m // 2, m)
    s = np.polyfit(q[sel, 0], q[sel, 2], 1)[0]; k0 = int(0.8 * m); J = q[k0:, 2].mean(); E = q[k0:, 6].mean()
    surv = (not stop) and abs(J - n) < 0.005 and abs(s) < 0.1 * abs(s0)
    cls = 'PRZETRWA' if surv else ('PEKA' if (stop or abs(J - n) >= 0.005) else 'POSREDNIE')
    r.update(ok=True, s=s, J=J, E=E, cls=cls)
    if q.shape[1] >= 15:
        dt = np.diff(q[:, 1]); rr = 0.5 * (q[1:, 9] + q[1:, 10]); plt_ = 2 * dt / (3 * c ** 3 * mu ** 2 * rr ** 4)
        wg = q[1:, 13]; wrr = q[1:, 14]; dE = np.diff(q[:, 6]); kk = max(int(0.8 * (m - 1)), 0)
        Erad = -wrr[kk:].mean(); PLT = plt_[kk:].mean(); Wg = wg[kk:].mean(); full = (wg[kk:] - dE[kk:]).mean()
        r.update(Erad=Erad, PLT=PLT, Wg=Wg, full=full)
    return r
R = {os.path.basename(p)[:-4]: summ(p) for p in sorted(glob.glob(os.path.join(O, '*.txt'))) if not os.path.basename(p).startswith('C')}
def row(r):
    s = f"  {r['name']:11s} {r['sp']:2s} {r['ch']:9s} n {r['n']} c {r['c']:8.3f} obiegi {r['orb']:6.0f} cykli {r['cyc']:5d}"
    if not r['ok']: return s + f"  -- {r['cls']}{' STOP fall' if r['fall'] else (' STOP' if r['stop'] else '')}"
    s += f"  J_k {r['J']:.6f} s {r['s']:+.3e} s/s0 {r['s'] / r['s0']:+.4f}  {r['cls']}{' STOP fall' if r['fall'] else (' STOP' if r['stop'] else '')}"
    if 'Erad' in r: s += f"  E_rad/(P_L T) {r['Erad'] / r['PLT']:.4f}  W_g/E_rad {r['Wg'] / r['Erad']:.4f}  [pelny bilans (W_g - dE)/(P_L T) {r['full'] / r['PLT']:.4f}]"
    return s
pr("=== wszystkie biegi 875 ==="); [pr(row(R[k])) for k in sorted(R)]
cc = [open(os.path.join(O, f'C{k}_871.txt'), 'rb').read() == open(os.path.join(O, f'C{k}_875.txt'), 'rb').read() for k in (1, 2)]
pr(f"H875c (LL=0 = local871 bajt w bajt: {cc}): {'TAK' if all(cc) else 'NIE'}")
K = [r for r in R.values() if r['name'].startswith('K_')]; okr = True
for sp in ('H', 'Ps'):
    sub = [r for r in K if r['sp'] == sp and r['ok']]; g = [0.8 <= r['s'] / r['s0'] <= 1.2 for r in sub]
    pr(f"  K {sp}: s/s0 {[(r['name'], round(r['s'] / r['s0'], 4)) for r in sub]} -> {sum(g)}/{len(sub)}; pominiete (< 10 cykli) {[r['name'] for r in K if r['sp'] == sp and not r['ok']]}")
    okr &= all(g) and len(sub) > 0
pr(f"H875r (kappa 0: s/s0 w [0,8; 1,2]): {'TAK' if okr else 'NIE'}")
F = {k: R[k] for k in R if k.startswith('F_')}
f137 = [F[f'F_{s}{n}_c137'] for s in 'HP' for n in (1, 2)]
pr(f"H875b (GLOWNA; c 137 PRZETRWA 4/4 -> {[(r['name'], r['cls']) for r in f137]}): {'TAK' if all(r['cls'] == 'PRZETRWA' for r in f137) else 'NIE'}")
S = [R[k] for k in R if k.startswith('S_')]; okm = True
for ch in ('F=0', 'F=1', 'p-Ps', 'o-Ps'):
    sub = [r for r in S if r['ch'] == ch]; g = sum(r['cls'] == 'PRZETRWA' for r in sub); pr(f"  {ch}: PRZETRWA {g}/{len(sub)} {[(r['name'], r['cls']) for r in sub]}"); okm &= (g == 2 and len(sub) == 2)
pr(f"H875m (kazdy kanal 2/2): {'TAK' if okm else 'NIE'}")
f10 = [F[f'F_{s}{n}_c10'] for s in 'HP' for n in (1, 2)]
pr(f"H875k (c 10 PEKA 4/4 -> {[(r['name'], r['cls']) for r in f10]}): {'TAK' if all(r['cls'] == 'PEKA' for r in f10) else 'NIE'}")
t1 = [F[f'F_{s}1_c30'] for s in 'HP']; t2 = [F[f'F_{s}2_c30'] for s in 'HP']
okt = all(r['cls'] == 'PEKA' for r in t1) and all(r['cls'] == 'PRZETRWA' for r in t2)
pr(f"H875t (c 30: n = 1 PEKA {[r['cls'] for r in t1]}, n = 2 PRZETRWA {[r['cls'] for r in t2]}): {'TAK' if okt else 'NIE'}")
sv = [r for r in list(F.values()) + S if r['c'] > 100 and r['cls'] == 'PRZETRWA']; okz = len(sv) > 0
for r in sv:
    a = r['Erad'] / r['PLT']; b = abs(r['Wg'] - r['Erad']) / r['Erad']; g = 0.8 <= a <= 1.2 and b < 0.2; okz &= g
    pr(f"  {r['name']:11s} {r['sp']:2s} {r['ch']:9s}: E_rad/(P_L T) {a:.4f}, |W_g - E_rad|/E_rad {b:.4f} -> {'T' if g else 'N'}; pelny bilans/(P_L T) {r['full'] / r['PLT']:.4f}")
pr(f"H875z (stan zablokowany promieniuje pelna moca, pole doplaca; {len(sv)} biegow): {'TAK' if okz else 'NIE'}")
meas = {'H': (2466061413.187e6, 'Parthey 2011'), 'Ps': (1233607216.4e6, 'Fee 1993')}; oke = all(r['cls'] == 'PRZETRWA' for r in f137)
if oke:
    for sp, key in (('H', 'H'), ('Ps', 'P')):
        r1, r2 = F[f'F_{key}1_c137'], F[f'F_{key}2_c137']; fm, cit = meas[sp]; f = (r2['E'] - r1['E']) * EH
        fj = EH * r1['mu'] * 0.5 * (1 / r1['J'] ** 2 - 1 / r2['J'] ** 2)
        pr(f"  {sp}: E1 {r1['E']:.9f} E2 {r2['E']:.9f} -> {f / 1e6:,.3f} MHz = {C0 / f * 1e9:.4f} nm = {f * HEV:.6f} eV; pomiar {fm / 1e6:,.3f} MHz ({cit}); "
           f"odch. {f / fm - 1:+.3e} = {(f - fm) * HEV * 1e6:+.1f} ueV  [z J: {fj / 1e6:,.3f} MHz, {fj / fm - 1:+.3e}]")
        oke &= abs(f / fm - 1) < 1e-4
pr(f"H875e (1S-2S z energii LOCK przy c 137, |odch.| < 1e-4): {'TAK' if oke else 'NIE (lub nie dotyczy)'}")

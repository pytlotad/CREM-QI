# Audit 871: local-field prototype -- verdicts H871c, b, k, l, e, r, d, m -> out/ana871.txt.  Written after prereg871.txt, before any result
# with kappa != 0 (only the kappa = 0 smoke tests were seen).
import os, glob, math, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
TWO_PI = 2 * math.pi
def load(name):
    p = os.path.join(O, name + '.txt'); hdr = {}; rows = []; stop = False
    for ln in open(p):
        if ln.startswith('# 871'):
            t = ln.split()[2:]; hdr = {t[i]: float(t[i + 1]) for i in range(0, len(t) - 1, 2)}
        elif ln.startswith('# STOP'): stop = True
        elif not ln.startswith('#') and ln.strip(): rows.append([float(x) for x in ln.split()])
    q = np.array(rows) if rows else np.zeros((0, 13)); return hdr, q, stop
def summ(name, kind='int'):
    hdr, q, stop = load(name); n = len(q)
    if n < 20: return dict(name=name, ok=False, n=n, stop=stop, fin=False, hdr=hdr)
    J = q[:, 2]; fin = bool(np.all(np.isfinite(q))) and not stop
    end = J[int(0.8 * n):].mean(); drift = abs(J[int(0.9 * n):].mean() - J[int(0.7 * n):int(0.8 * n)].mean())
    ms1, ms2 = hdr.get('MS1', 0.0), hdr.get('MS2', 0.0); P1 = q[int(0.8 * n):, 7].mean(); P2 = q[int(0.8 * n):, 8].mean()
    shift = (ms1 * P1 + ms2 * P2) / TWO_PI
    if kind == 'int': tgt = round(end - shift) + shift
    elif kind == '2/3': tgt = round(end * 1.5) / 1.5
    lock = drift < 0.002 and abs(end - tgt) < 0.005 and fin
    return dict(name=name, ok=True, n=n, J=end, drift=drift, tgt=tgt, dev=end - tgt, lock=lock, conv=drift < 0.002 and fin, fin=fin, stop=stop,
                L=q[int(0.8 * n):, 4].mean(), Lmin=q[:, 4].min(), rmin=q[:, 9].min(), E=q[int(0.8 * n):, 6].mean(), dth=q[int(0.8 * n):, 5].mean(),
                shift=shift, A0=q[0, 3], A1=q[-1, 3], E0=q[0, 6], E1=q[-1, 6], hdr=hdr, P1=P1, P2=P2)
def table(rows, title, kind='int'):
    pr(f"=== {title} ===");  pr("  bieg               cykle  J_koniec      dryf      cel        odch.      LOCK  L         min L     min r     E")
    for r in rows:
        if not r['ok']: pr(f"  {r['name']:18s} {r['n']:5d}  (za malo cykli{', STOP' if r['stop'] else ''})"); continue
        pr(f"  {r['name']:18s} {r['n']:5d}  {r['J']:.6f}  {r['drift']:.2e}  {r['tgt']:.6f}  {r['dev']:+.2e}  {'T' if r['lock'] else 'N'}     {r['L']:+.5f}  {r['Lmin']:+.5f}  {r['rmin']:.4f}  {r['E']:+.6e}")
def names(prefix): return sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(O, prefix + '*.txt')))
# H871c
okc = True
for nm in ('K_c', 'K_e'):
    hdr, q, stop = load(nm); da = abs(q[-1, 3] / q[0, 3] - 1); de = abs(q[-1, 6] / q[0, 6] - 1); okc &= da < 1e-6 and de < 1e-6 and not stop
    pr(f"  {nm}: cykli {len(q)}, |A_koniec/A_start - 1| {da:.1e}, |E_koniec/E_start - 1| {de:.1e}")
pr(f"H871c (kappa 0: A i E stale do 1e-6): {'TAK' if okc else 'NIE'}")
# H871b
def false_lock(r):
    if not r['conv']: return False
    J = r['J']; k = round(2 * J)
    if k % 2 == 1 and abs(J - k / 2) < 0.005: return True
    k3 = round(1.5 * J); v = k3 / 1.5
    return abs(J - v) < 0.005 and abs(v - round(v)) > 0.01
okb = True; G = {}
for pre, nm in (('GH_', 'wodor'), ('GP_', 'pozytonium')):
    R = [summ(x) for x in names(pre)]; G[nm] = R; table(R, f"(G{pre[1]}) {nm}, S3, pole lokalne, nierelatywistycznie")
    nl = sum(r.get('lock', False) and abs(r['tgt'] - round(r['tgt'])) < 1e-9 for r in R); nf = sum(r['ok'] and false_lock(r) for r in R)
    pr(f"  {nm}: LOCK przy calkowitym N {nl}/{len(R)}, falszywe zbieznosci (polcalkowite / 2k/3) {nf}")
    okb &= nl >= 15 and nf == 0
pr(f"H871b (GLOWNA; >= 15/18 LOCK calkowite i 0 falszywych, osobno wodor i pozytonium): {'TAK' if okb else 'NIE'}")
# H871k
R = [summ(x, '2/3') for x in names('GL_')]; table(R, '(GL) wodor, zegar S1 (lagranzjan), cel 2k/3', '2/3')
conv = [r for r in R if r['ok'] and r['conv']]; at = [r for r in conv if abs(r['dev']) < 0.005]; nonint = [r for r in at if abs(r['tgt'] - round(r['tgt'])) > 0.01]
pr(f"  zbiezne {len(conv)}/{len(R)}, przy 2k/3 {len(at)}, w tym niecalkowite {len(nonint)}: {[round(r['J'], 4) for r in nonint]}")
okk = len(conv) > 0 and len(at) >= 0.8 * len(conv) and len(nonint) >= 2
pr(f"H871k (>= 80 % zbieznych przy 2k/3, >= 2 niecalkowite): {'TAK' if okk else 'NIE'}")
# H871l
R = [summ(x) for x in names('GZ')]; table(R, '(GZ) e 0,95 (orbity prawie radialne), S3')
okl = all(r['ok'] and r['fin'] and r['rmin'] > 0 and r['Lmin'] > 0 for r in R) and sum(r.get('lock', False) for r in R) >= 6
pr(f"  LOCK {sum(r.get('lock', False) for r in R)}/{len(R)}; min L {min(r['Lmin'] for r in R if r['ok']):+.5f}; min r {min(r['rmin'] for r in R if r['ok']):.4f}")
pr(f"H871l (skonczone, min r > 0, min L > 0 w 8/8 i LOCK >= 6/8): {'TAK' if okl else 'NIE'}")
# H871e
RINF = 10973731.568157 * 299792458.0; MEMP = 5.446170214889e-4; C0 = 299792458.0; HEV = 4.135667696e-15
meas = {'wodor': (2466061413.187e6, 1 / (1 + MEMP), 'Parthey 2011', 1 / (1 + MEMP)), 'pozytonium': (1233607216.4e6, 0.5, 'Fee 1993', 0.5)}
oke = True
pr("(H871e) 1S-2S z LOCK N = 1 i N = 2 (poziomy grube bez spinu: F = 0 / F = 1 i p-Ps / o-Ps jednakowe):")
for nm, (fm, red, cit, frac1) in meas.items():
    L1 = [r['J'] for r in G[nm] if r.get('lock') and round(r['tgt']) == 1]; L2 = [r['J'] for r in G[nm] if r.get('lock') and round(r['tgt']) == 2]
    if not L1 or not L2: pr(f"  {nm}: brak LOCK przy N = 1 ({len(L1)}) albo N = 2 ({len(L2)})"); oke = False; continue
    J1, J2 = np.mean(L1), np.mean(L2); f = RINF * red * (1 / J1 ** 2 - 1 / J2 ** 2)
    # per-particle closure (electron alone): oint p_2 dx_2 = (M1/MT) oint p dr = 2 pi n -> N = n / frac(M1/MT)
    fpp = RINF * red * (frac1 ** 2) * (1 - 1 / 4)
    pr(f"  {nm}: J1 {J1:.6f} ({len(L1)}), J2 {J2:.6f} ({len(L2)}) -> {f / 1e6:.3f} MHz = {f * HEV:.4f} eV, {C0 / f * 1e9:.3f} nm; pomiar {fm / 1e6:.3f} MHz ({cit}); odch. {f / fm - 1:+.3e}")
    pr(f"     domkniecie osobno dla czastki 2: {fpp / 1e6:.3f} MHz, odch. {fpp / fm - 1:+.3e}")
    oke &= abs(f / fm - 1) < 0.01 and abs(fpp / fm - 1) > (1e-4 if nm == 'wodor' else 0.10)
pr(f"H871e (LOCK 1S-2S w 1 %; domkniecie osobno: wodor > 1e-4, Ps > 10 %): {'TAK' if oke else 'NIE'}")
# H871r, H871d
C = 10.0; A = 1 / C
def somm(Jr, L, mu): return mu * C * C * ((1 + A * A / (Jr + math.sqrt(max(L * L - A * A, 0))) ** 2) ** -0.5 - 1)
RH = {x[3:]: summ(x) for x in names('RH_')}; RP = {x[3:]: summ(x) for x in names('RP_')}; RN = {x[3:]: summ(x) for x in names('RN_')}
table(RH.values(), '(RH) wodor, LW c = 10, opozniony zapis'); table(RP.values(), '(RP) pozytonium, LW c = 10, opozniony zapis (raport)'); table(RN.values(), '(RN) wodor, LW bez opoznionego zapisu')
pr("  RH: energia wobec Sommerfelda klasycznie E(J_r = J - |L|, |L|) (mu = m_p/(m_p+1)):")
for k, r in RH.items():
    if r['ok']: mu = 1 / (1 + MEMP); Es = somm(r['J'] - abs(r['L']), abs(r['L']), mu); Eb = -mu / (2 * r['J'] ** 2)
    if r['ok']: pr(f"    {k}: E {r['E']:+.6e}, E_Somm {Es:+.6e}, E_Bohr(J) {Eb:+.6e}, (E - E_S)/(E_S - E_B) {(r['E'] - Es) / (Es - Eb):+.3f}; L {r['L']:+.4f} (calkowite? {abs(abs(r['L']) - round(abs(r['L']))) < 0.01})")
okr = sum(r.get('lock', False) and abs(r['tgt'] - round(r['tgt'])) < 1e-9 for r in RH.values()) >= 5
pr(f"H871r (RH LOCK J_r + L calkowite >= 5/6): {'TAK' if okr else 'NIE'}")
pairs = [(k, RH[k], RN[k]) for k in RH if k in RN and RH[k].get('lock') and RN[k].get('lock')]
for k, a, b in pairs: pr(f"    {k}: J(RH) - J(RN) {a['J'] - b['J']:+.2e}")
okd = len(pairs) > 0 and all(abs(a['J'] - b['J']) < 1e-3 for _, a, b in pairs)
pr(f"H871d (opozniony zapis bez wplywu, |dJ| < 1e-3, par {len(pairs)}): {'TAK' if okd else 'NIE'}")
# H871m
okm = True; CH = {'SH0': 'wodor F = 0', 'SH1': 'wodor F = 1', 'SP0': 'p-Ps', 'SP1': 'o-Ps'}; SM = {}
for pre, nm in CH.items():
    R = [summ(x) for x in names(pre + '_')]; SM[pre] = {x['name'][4:]: x for x in R}; table(R, f"(S) {nm}: cel n + sum m Phi / 2 pi")
    good = 0
    for r in R:
        if not r['ok']: continue
        tol = 2e-4 if abs(r['shift']) < 2e-3 else 0.1 * abs(r['shift']); g = r['lock'] and abs(r['dev']) <= tol; good += g
        pr(f"    {r['name']}: przesuniecie {r['shift']:+.3e} (Phi1 {r['P1']:+.3e}, Phi2 {r['P2']:+.3e}), odch. {r['dev']:+.2e}, tolerancja {tol:.1e} -> {'T' if g else 'N'}")
    pr(f"  {nm}: {good}/{len(R)}"); okm &= good >= 3
pr(f"H871m (kazdy kanal >= 3/4 przy celu z faza spinu): {'TAK' if okm else 'NIE'}")
pr("  roznice kanalow (te same starty; jednostki prototypu, c = 10 -- porownanie z pomiarem NIE DOTYCZY):")
for a, b, nm in (('SH1', 'SH0', 'wodor F=1 - F=0'), ('SP1', 'SP0', 'o-Ps - p-Ps')):
    for k in SM[a]:
        if k in SM[b] and SM[a][k]['ok'] and SM[b][k]['ok']:
            pr(f"    {nm} {k}: dJ {SM[a][k]['J'] - SM[b][k]['J']:+.3e} (przewid. {SM[a][k]['shift'] - SM[b][k]['shift']:+.3e}), dE {SM[a][k]['E'] - SM[b][k]['E']:+.3e}")

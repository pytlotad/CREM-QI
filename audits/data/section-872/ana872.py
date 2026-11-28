# Audit 872: local-field prototype, starts in the capture window, TEND 6e4 -- verdicts H872c, b, p, k, l, e, r, s, m, f -> out/ana872.txt.
# Written after prereg872.txt and the launch of the runs, before any result was looked at.
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
def summ(name, kind='n'):
    hdr, q, stop = load(name); n = len(q)
    if n < 20: return dict(name=name, ok=False, n=n, stop=stop, fin=False, hdr=hdr, lock=False, conv=False)
    J = q[:, 2]; fin = bool(np.all(np.isfinite(q))) and not stop
    end = J[int(0.8 * n):].mean(); drift = abs(J[int(0.9 * n):].mean() - J[int(0.7 * n):int(0.8 * n)].mean())
    ms1, ms2 = hdr.get('MS1', 0.0), hdr.get('MS2', 0.0); P1 = q[int(0.8 * n):, 7].mean(); P2 = q[int(0.8 * n):, 8].mean()
    shift = (ms1 * P1 + ms2 * P2) / TWO_PI
    if kind == 'n': tgt = round(hdr['N0']) + shift
    else: tgt = round(end * 1.5) / 1.5
    conv = drift < 0.002 and fin
    lock = conv and abs(end - tgt) < 0.005
    return dict(name=name, ok=True, n=n, J=end, drift=drift, tgt=tgt, dev=end - tgt, lock=lock, conv=conv, fin=fin, stop=stop,
                L=q[int(0.8 * n):, 4].mean(), Lmin=q[:, 4].min(), rmin=q[:, 9].min(), E=q[int(0.8 * n):, 6].mean(),
                shift=shift, A0=q[0, 3], A1=q[-1, 3], E0=q[0, 6], E1=q[-1, 6], hdr=hdr, P1=P1, P2=P2, N0=hdr['N0'], e=hdr['E0'])
def table(rows, title):
    pr(f"=== {title} ===");  pr("  bieg               cykle  J_koniec      dryf      cel        odch.      LOCK  L         min L     min r     E")
    for r in rows:
        if not r['ok']: pr(f"  {r['name']:18s} {r['n']:5d}  (za malo cykli{', STOP' if r['stop'] else ''})"); continue
        pr(f"  {r['name']:18s} {r['n']:5d}  {r['J']:.7f}  {r['drift']:.2e}  {r['tgt']:.6f}  {r['dev']:+.2e}  {'T' if r['lock'] else 'N'}     {r['L']:+.5f}  {r['Lmin']:+.5f}  {r['rmin']:.4f}  {r['E']:+.7e}{'  STOP' if r['stop'] else ''}")
def names(prefix): return sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(O, prefix + '*.txt')))
# H872c
hdr, q, stop = load('K_c'); da = abs(q[-1, 3] / q[0, 3] - 1); de = abs(q[-1, 6] / q[0, 6] - 1); okc = da < 1e-6 and de < 1e-6 and not stop
pr(f"  K_c: cykli {len(q)}, t {q[-1, 1]:.0f}, |A_koniec/A_start - 1| {da:.1e}, |E_koniec/E_start - 1| {de:.1e}")
pr(f"H872c (kappa 0: A i E stale do 1e-6): {'TAK' if okc else 'NIE'}")
# H872b, H872p
okb = okp = True; G = {}
for pre, nm in (('GH_', 'wodor'), ('GP_', 'pozytonium')):
    R = [summ(x) for x in names(pre)]; G[nm] = R; table(R, f"(G{pre[1]}) {nm}, S3, starty w oknie, TEND 6e4")
    nl = sum(r['lock'] for r in R); n2 = [r for r in R if r['ok'] and abs(abs(r['N0'] - round(r['N0'])) - 0.02) < 1e-9]
    nl2 = sum(r['lock'] for r in n2); nf = sum(r['conv'] and not r['lock'] for r in R)
    for n in (1, 2, 3):
        sub = [r for r in R if r['ok'] and round(r['N0']) == n]
        pr(f"  {nm} n = {n}: LOCK {sum(r['lock'] for r in sub)}/{len(sub)}")
    devs = [abs(r['dev']) for r in R if r['lock']]
    med = float(np.median(devs)) if devs else float('nan')
    pr(f"  {nm}: LOCK {nl}/{len(R)}, przy |d| = 0,02 {nl2}/{len(n2)}, falszywe {nf}; mediana |J - n| po LOCK {med:.2e} (zakres {min(devs) if devs else 0:.1e} - {max(devs) if devs else 0:.1e})")
    okb &= nl >= 20 and nl2 >= 11 and nf == 0; okp &= bool(devs) and med < 1e-4
pr(f"H872b (GLOWNA; >= 20/24, |d| 0,02 >= 11/12, 0 falszywych; osobno wodor i pozytonium): {'TAK' if okb else 'NIE'}")
pr(f"H872p (mediana |J - n| po LOCK < 1e-4, osobno): {'TAK' if okp else 'NIE'}")
# H872k
R = [summ(x, '2/3') for x in names('GL_')]; table(R, '(GL) wodor, zegar S1 (lagranzjan); cel = najblizsza wielokrotnosc 2/3')
okk = True
for r in R:
    N0 = r['N0']; mid = abs(N0 - round(N0)) < 0.05 and round(N0) in (1, 3)
    if mid:
        bad = r['ok'] and r['conv'] and abs(r['J'] - round(r['J'])) < 0.005 and round(r['J']) % 2 == 1
        pr(f"  {r['name']}: start przy {round(N0)} (polowa miedzy 2k/3) -> J {r.get('J', float('nan')):.5f}, LOCK przy calkowitym? {'TAK (zle)' if bad else 'nie'}"); okk &= not bad
    else:
        want = round(N0 * 1.5) / 1.5; g = r['ok'] and r['lock'] and abs(r['tgt'] - want) < 1e-9
        pr(f"  {r['name']}: cel {want:.4f} -> LOCK {'T' if g else 'N'} (J {r.get('J', float('nan')):.6f}, odch. {r.get('dev', float('nan')):+.2e})"); okk &= g
pr(f"H872k (1 i 3 bez LOCK calkowitego; 0,69 / 1,30 LOCK 2/3 / 4/3): {'TAK' if okk else 'NIE'}")
# H872l
R = [summ(x) for x in names('GZ')]; table(R, '(GZ) e 0,95 (orbity prawie radialne), S3')
okl = all(r['ok'] and r['fin'] and r['rmin'] > 0 and r['Lmin'] > 0 for r in R) and sum(r['lock'] for r in R) >= 6
pr(f"  LOCK {sum(r['lock'] for r in R)}/{len(R)}; min L {min(r['Lmin'] for r in R if r['ok']):+.5f}; min r {min(r['rmin'] for r in R if r['ok']):.4f}")
pr(f"H872l (skonczone, min r > 0, min L > 0 w 8/8 i LOCK >= 6/8): {'TAK' if okl else 'NIE'}")
# H872e
RINF = 10973731.568157 * 299792458.0; MEMP = 5.446170214889e-4; C0 = 299792458.0; HEV = 4.135667696e-15
meas = {'wodor': (2466061413.187e6, 1 / (1 + MEMP), 'Parthey 2011'), 'pozytonium': (1233607216.4e6, 0.5, 'Fee 1993')}
oke = True
pr("(H872e) 1S-2S z LOCK n = 1 i n = 2 (poziomy grube bez spinu: F = 0 / F = 1 i p-Ps / o-Ps jednakowe):")
for nm, (fm, red, cit) in meas.items():
    L1 = [r['J'] for r in G[nm] if r['lock'] and round(r['N0']) == 1]; L2 = [r['J'] for r in G[nm] if r['lock'] and round(r['N0']) == 2]
    if not L1 or not L2: pr(f"  {nm}: brak LOCK przy n = 1 ({len(L1)}) albo n = 2 ({len(L2)})"); oke = False; continue
    J1, J2 = np.mean(L1), np.mean(L2); f = RINF * red * (1 / J1 ** 2 - 1 / J2 ** 2); fb = RINF * red * 0.75
    pr(f"  {nm}: J1 {J1:.7f} (sd {np.std(L1):.1e}, {len(L1)}), J2 {J2:.7f} (sd {np.std(L2):.1e}, {len(L2)})")
    pr(f"     -> {f / 1e6:,.3f} MHz = {C0 / f * 1e9:.4f} nm = {f * HEV:.6f} eV; pomiar {fm / 1e6:,.3f} MHz ({cit}); odch. {f / fm - 1:+.3e} = {(f - fm) / 1e6:+,.1f} MHz = {(f - fm) * HEV * 1e6:+.1f} ueV")
    pr(f"     Bohr (J calkowite) {fb / 1e6:,.3f} MHz, odch. od pomiaru {fb / fm - 1:+.3e}; prototyp wobec Bohra {f / fb - 1:+.3e}")
    oke &= abs(f / fm - 1) < 1e-3
pr(f"H872e (1S-2S z LOCK w 1e-3 od pomiaru, osobno): {'TAK' if oke else 'NIE'}")
# H872r
C = 10.0; A = 1 / C
def somm(Jr, L, mu): return mu * C * C * ((1 + A * A / (Jr + math.sqrt(max(L * L - A * A, 0))) ** 2) ** -0.5 - 1)
RH = [summ(x) for x in names('RH_')]; table(RH, '(RH) wodor, LW c = 10, opozniony zapis')
pr("  RH: energia wobec Sommerfelda klasycznie E(J_r = J - |L|, |L|):")
for r in RH:
    if r['ok']:
        mu = 1 / (1 + MEMP); Es = somm(r['J'] - abs(r['L']), abs(r['L']), mu); Eb = -mu / (2 * r['J'] ** 2)
        pr(f"    {r['name']}: E {r['E']:+.7e}, E_Somm {Es:+.7e}, E_Bohr(J) {Eb:+.7e}, (E - E_S)/(E_S - E_B) {(r['E'] - Es) / (Es - Eb):+.3f}; L {r['L']:+.5f} (calkowite w 0,01? {abs(abs(r['L']) - round(abs(r['L']))) < 0.01})")
okr = sum(r['lock'] for r in RH) >= 5
pr(f"H872r (RH LOCK przy n >= 5/6): {'TAK' if okr else 'NIE'}")
# H872s
RP = [summ(x) for x in names('RP_')]; table(RP, '(RP) pozytonium, LW c = 30'); oks = sum(r['lock'] for r in RP) >= 3
pr(f"H872s (RP LOCK >= 3/4): {'TAK' if oks else 'NIE'}")
RQ = [summ(x) for x in names('RQ_')]; table(RQ, '(RQ) pozytonium, LW c = 10 (raport)')
# H872m, H872f
okm = True; CH = {'SH0': 'wodor F = 0 (c 10)', 'SH1': 'wodor F = 1 (c 10)', 'SP0': 'p-Ps (c 30)', 'SP1': 'o-Ps (c 30)'}; SM = {}
for pre, nm in CH.items():
    R = [summ(x) for x in names(pre + '_')]; SM[pre] = {x['name'][4:]: x for x in R}; table(R, f"(S) {nm}: cel n + sum m Phi / 2 pi")
    good = 0
    for r in R:
        if not r['ok']: continue
        tol = 2e-4 if abs(r['shift']) < 2e-3 else 0.1 * abs(r['shift']); g = r['lock'] and abs(r['dev']) <= tol; good += g
        pr(f"    {r['name']}: przesuniecie {r['shift']:+.3e} (Phi1 {r['P1']:+.3e}, Phi2 {r['P2']:+.3e}), odch. {r['dev']:+.2e}, tolerancja {tol:.1e} -> {'T' if g else 'N'}")
    pr(f"  {nm}: {good}/{len(R)}"); okm &= good >= 3
pr(f"H872m (kazdy kanal >= 3/4 przy celu z faza momentow): {'TAK' if okm else 'NIE'}")
pr("  roznice kanalow (te same starty, oba LOCK; jednostki prototypu -- porownanie z pomiarem NIE DOTYCZY):")
okf = True; nf = 0
for a, b, nm in (('SH1', 'SH0', 'wodor F=1 - F=0'), ('SP1', 'SP0', 'o-Ps - p-Ps')):
    for k in SM[a]:
        if k in SM[b] and SM[a][k]['lock'] and SM[b][k]['lock']:
            dj = SM[a][k]['J'] - SM[b][k]['J']; dt = SM[a][k]['shift'] - SM[b][k]['shift']
            if abs(dt) < 2e-6: pr(f"    {nm} {k}: dJ {dj:+.3e}, roznica celow {dt:+.3e} < 2e-6 -> NIE DOTYCZY"); continue
            g = np.sign(dj) == np.sign(dt) and abs(dj / dt - 1) <= 0.3; okf &= g; nf += 1
            pr(f"    {nm} {k}: dJ {dj:+.3e}, roznica celow {dt:+.3e}, stosunek {dj / dt:+.3f}, dE {SM[a][k]['E'] - SM[b][k]['E']:+.3e} -> {'T' if g else 'N'}")
pr(f"H872f (znak i 30 %, porownan {nf}): {'TAK' if okf and nf > 0 else 'NIE'}")

# Audit 786: where the clock is read -- proper time on the particle vs spatial phase in the centre-of-mass frame vs lab frame (boost) -> out/ana786.txt
import os, glob, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
pr = lambda *a: print(*a, flush=True)
h = 6.62607015e-34; ME = 9.1093837139e-31; MP = 1.67262192595e-27; ALPHA = 7.2973525643e-3
MASS = {'H': (MP, ME), 'aH': (MP, ME), 'pPs': (ME, ME), 'oPs': (ME, ME)}
def parse(fn):
    start, per = None, []
    for l in open(fn):
        if l.startswith('PAIR CLOCK701 start'):
            p = l.split(); start = dict(nlock=float(p[p.index('nlock') + 1]), binding=float(p[p.index('binding') + 1]))
        elif l.startswith('PAIR CLOCK701 k '):
            p = l.split()[2:]; per.append({p[i]: float(p[i + 1]) for i in range(0, len(p) - 1, 2)})
    assert start and len(per) == 20, (fn, len(per))
    return start, per
SYS = ['H', 'aH', 'pPs', 'oPs']; SEEDS = [0, 1]
R = {(v, s, sd): parse(os.path.join(O, f"{v}_{s}_{sd}.err")) for v in ('S', 'BP', 'BA') for s in SYS for sd in SEEDS}
okA = okB = okC = okD = True
pr("=== S (bez boostu): odczyty zegara, mediany z 20 okresow ===")
pr("  uklad ziarno  J/h        J/h - n_lock   J/h - n_E      Tau/J        J1/J (oczek.)            Tau1/Tau")
for s in SYS:
    for sd in SEEDS:
        st, P = R[('S', s, sd)]; m1, m2 = MASS[s]; w = m2 / (m1 + m2)
        J = np.array([p['J'] for p in P]); tau = np.array([p['Tau1'] + p['Tau2'] for p in P])
        j1 = np.array([p['J1'] for p in P]) / J; t1 = np.array([p['Tau1'] for p in P]) / tau
        nE = np.sqrt(st['binding'] / -np.array([p['E'] for p in P]))
        a = np.median(tau / J); okA &= abs(a - 0.5) < 1e-4
        okB &= abs(np.median(j1) - w) < 1e-4 and abs(np.median(t1) - w) < 1e-4
        dl = np.max(np.abs(J / h - st['nlock'])); okD &= dl <= 1e-3
        pr(f"  {s:4s} {sd}     {np.median(J / h):.7f}  {np.median(J / h - st['nlock']):+.3e}     {np.median(J / h - nE):+.3e}     {a:.7f}    {np.median(j1):.7f} ({w:.7f})   {np.median(t1):.7f}")
pr(f"H786a (|Tau/J - 1/2| < 1e-4, wszystkie S): {'TAK' if okA else 'NIE'}")
pr(f"H786b (J1/J i Tau1/Tau = m2/M do 1e-4): {'TAK' if okB else 'NIE'}")
pr(f"H786d (|J/h - n_lock| <= 1e-3 w kazdym okresie S): {'TAK' if okD else 'NIE'}")
pr("=== BP, BA (boost 1e-3 c): (J_lab - PV) / J_S - 1, mediany ===")
for s in SYS:
    for sd in SEEDS:
        JS = np.median([p['J'] for p in R[('S', s, sd)][1]]); row = []
        for v in ('BP', 'BA'):
            P = R[(v, s, sd)][1]; d = np.median([(p['J'] - p['PV']) for p in P]) / JS - 1; okC &= abs(d) < 1e-4
            tl = np.median([(p['Tau1'] + p['Tau2']) / p['J'] for p in P])
            row.append(f"{v} {d:+.3e} (J_lab/h {np.median([p['J'] for p in P]) / h:.5f}, Tau/J_lab {tl:.6f})")
        pr(f"  {s:4s} {sd}   " + '   '.join(row))
pr(f"H786c (|(J_lab - PV)/J_S - 1| < 1e-4, BP i BA): {'TAK' if okC else 'NIE'}")
pr(f"  odniesienie: alpha^2 / 4 n^2 (n = 2) = {ALPHA ** 2 / 16:.2e}")

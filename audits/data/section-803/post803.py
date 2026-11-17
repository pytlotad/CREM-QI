# Audyt 803 POST HOC: miara na obieg odporna na falszywe apocentrum: srednia delta tylko z odcinkow |dt/mediana - 1| <= 0,1, / (2 pi T_n).
# (Wczesniejsze warianty tego skryptu, ujawnione w uzupelnieniu: sum(delta)/sum(t) -- miesza okres rzeczywisty; sum(delta)/round(sum dt/mediana) -- przestrzelil.)
import math, numpy as np
src = open('ana803.py').read(); head = src.split('W = {}; okV')[0]; ns = {'__file__': 'ana803.py'}; exec(head, ns)
sec, Ef = ns['sec'], ns['Ef']
def robust(p):
    T = []; L = []; S = []; per = None
    for ln in open(p):
        if ln.startswith('PAIR TORQUE799 mechB0 sum'): per = float(ln.split()[ln.split().index('period') + 1])
        if ln.startswith('PAIR TORQUE799 mechB0'):
            w = ln.split()
            if w[3] == 'apo': T.append(float(w[w.index('t') + 1])); L.append(ns['v3'](w, 'L')); S.append(ns['v3'](w, 's2'))
    d = [np.linalg.norm(S[j] - S[j - 1])/np.linalg.norm(np.cross(L[j - 1], S[j - 1])) for j in range(1, len(S))]
    dt = np.diff(T); ok = np.abs(dt/np.median(dt) - 1) <= 0.1
    if not ok.all(): print(f'  [{p.split("/")[-1]}: odrzucone {int((~ok).sum())} z {len(d)} odcinkow, dt/mediana ' + ', '.join('%.4f' % x for x in dt[~ok]/np.median(dt)) + ']')
    return float(np.mean(np.array(d)[ok]))/(2*math.pi*per)
TORI = [(4, x) for x in (1, 1.5, 2, 2.5, 3, 3.5, 4)]
for F, seeds in ((1, (4, 5, 6)), (0, (5,))):
    W = {}
    for s in seeds:
        for t in ('proton', 'antiproton'):
            for k in TORI:
                p = f'out/runs/M/F{F}_{t}_n4L{k[1]:g}_{s}.err'; W[(t, k, s)] = dict(obr=robust(p), S=sec(p)[0]['S'])
    pm = lambda k, s, x: 0.5*(W[('proton', k, s)][x] + W[('antiproton', k, s)][x])
    print(f'=== F = {F} (PHEN={2 if F else 1}), miara: srednia delta z czystych odcinkow / T_n ===')
    for l, st in ((1, '4p'), (2, '4d'), (3, '4f')):
        E = []
        for s in seeds:
            od = l*pm((4, l), s, 'obr') - (l + 1)*pm((4, l + 1), s, 'obr'); E.append(od/(pm((4, l + 0.5), s, 'obr')*(l + 0.5)**2/(l*(l + 1))) - 1)
        eS = np.mean([(l*pm((4, l), s, 'S') - (l + 1)*pm((4, l + 1), s, 'S'))/(pm((4, l + 0.5), s, 'S')*(l + 0.5)**2/(l*(l + 1))) - 1 for s in seeds])
        print(f'  {st}: E_mech {np.mean(E):+.3e} (ziarna {", ".join("%+.2e" % e for e in E)}) | E_S {eS:+.3e} | E_mech - E_S {np.mean(E) - eS:+.2e} | E_mech - E_f {np.mean(E) - Ef(l):+.2e}')
    for k in TORI:
        r = [pm(k, s, 'obr')/pm(k, s, 'S') - 1 for s in seeds]
        print(f'  n4L{k[1]:g}: obr/S - 1 (pary) ' + ', '.join('%+.2e' % x for x in r))

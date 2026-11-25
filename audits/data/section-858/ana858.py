# Audit 858: ORBIT858 trace of the transit run vs the orbit-averaged dipole equations.  Usage: LC_ALL=C python3 ana858.py -> out/ana858.txt.
import os
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
R = os.environ.get('R858', SP + '/r858'); pr = lambda *a: print(*a, flush=True)
rd = lambda v, t, x: open(f'{R}/reg/{v}/{t}/{x}', 'rb').read()
strip = lambda b: b'\n'.join(l for l in b.split(b'\n') if not l.startswith(b'ORBIT858'))
k1 = all(rd('old', t, x) == rd('new', t, x) for t in ('K_p2_n2', 'K_p1_n1') for x in ('out.txt', 'err.txt', 'rc.txt'))
pairs = (('K_p2_n2', 'P_p2_n2'), ('K_p1_n1', 'P_p1_n1'))
k2 = all(rd('old', a, 'out.txt') == rd('new', b, 'out.txt') and strip(rd('new', b, 'err.txt')) == rd('old', a, 'err.txt')
         and rd('old', a, 'rc.txt') == rd('new', b, 'rc.txt') for a, b in pairs)
pr(f"H858k (bez przelacznika new = old; z przelacznikiem out = old, err = old bez ORBIT858): {'TAK' if k1 and k2 else 'NIE'}  (k1 {k1}, k2 {k2})")
def parse(t):
    O, T = {}, {}
    for ln in rd('new', t, 'err.txt').decode(errors='replace').split('\n'):
        w = ln.split()
        if ln.startswith('ORBIT858'):
            d = {w[i]: w[i + 1] for i in range(1, len(w) - 1, 2)}; O.setdefault(d['seed'], []).append(d)
        elif ln.startswith('TRANSIT855'):
            d = {w[i]: w[i + 1] for i in range(1, len(w) - 1, 2)}; T[d['seed']] = d
    return O, T
f = lambda d, k: float(d[k])
okZ, okC, okE, okM = True, True, True, True
for p, lab in (('1', 'p-Ps (PHEN=1)'), ('2', 'o-Ps (PHEN=2)')):
    for kind, n in (('P', '1'), ('P', '2'), ('C', '2')):
        t = f'{kind}_p{p}_n{n}'; O, T = parse(t)
        pr(f'  {kind} {lab} n = {n}: ziaren {len(O)}, TRANSIT855 {len(T)}, rc {rd("new", t, "rc.txt").decode().strip()}')
        for s, L in O.items():
            last = L[-1]; tl = f(T[s], 'tlab') if s in T else float('nan')
            z = abs(f(last, 'tlab') / tl - 1) if s in T else float('nan')
            if not z <= 1e-6: okZ = False
            dt = [f(L[0], 't')] + [f(L[i], 't') - f(L[i - 1], 't') for i in range(1, len(L))]
            wl = sum(f(d, 'loss') * x for d, x in zip(L, dt)) / sum(dt)
            pers = [f(d, 'per') for d in L]
            dev = [f(d, 'e1') / f(d, 'einv') - 1 for d in L if f(d, 'einv') >= 0.05]
            pr(f'    seed {s}: ck {len(L)}, tlab/TRANSIT855 - 1 = {z:.2e}, loss (wazona) {wl:.5f} [{min(f(d, "loss") for d in L):.5f}, {max(f(d, "loss") for d in L):.5f}],'
               f' per [{min(pers):.6f}, {max(pers):.6f}], max |e/einv - 1| {max((abs(x) for x in dev), default=0):.4f}, t/tavg (koniec) {f(last, "t") / f(last, "tavg"):.4f}')
            for d in L[:3] + L[-1:]:
                pr(f"      ck {d['ck']}: a0 {f(d, 'a0'):.4f} -> a1 {f(d, 'a1'):.4f}, e {f(d, 'e0'):.4f} -> {f(d, 'e1'):.4f} (einv {f(d, 'einv'):.4f}), skip {d['skip']}, loss {f(d, 'loss'):.5f}, t/tavg {f(d, 't') / f(d, 'tavg'):.4f}")
            if kind == 'C' and not (0.99 <= wl <= 1.01 and 0.999 <= min(pers) and max(pers) <= 1.001): okC = False
            if kind == 'P' and n == '2':
                if not all(abs(x) <= 0.02 for x in dev): okE = False
                if not 1.06 <= wl <= 1.12: okM = False
pr(f"H858z (ostatni ORBIT858 tlab = TRANSIT855 tlab do 1e-6): {'TAK' if okZ else 'NIE'}")
pr(f"H858c (start kolowy n = 2: loss wazona w [0,99; 1,01], per w [0,999; 1,001]): {'TAK' if okC else 'NIE'}")
pr(f"H858e (n = 2 domyslnie: |e/einv - 1| <= 0,02 przy einv >= 0,05): {'TAK' if okE else 'NIE'}")
pr(f"H858m (GLOWNA; n = 2 domyslnie: loss wazona w [1,06; 1,12], p-Ps i o-Ps): {'TAK' if okM else 'NIE'}")

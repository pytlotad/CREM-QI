# Audit 819 (repeat of 818, new seeds): default Plummer floor 0.05 -> 0.01 r*.  Written during the runs, before reading them.  Usage: python3 ana819.py -> out/ana819.txt
import os, sys, re, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); sys.path[:0] = [D, os.path.join(D, '..', 'section-815'), os.path.join(D, '..', 'section-814')]
import pred815 as P15
from pred814 import parse_start
O = os.path.join(D, 'out'); pr = lambda *a: print(*a, flush=True)
HEV, C = 4.135667696e-15, 299792458.0
units = lambda f: f'{f/1e6:.4f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
# ---------------- H819o, H819h (byte-for-byte)
reg = dict((tuple(l.split()[2:4]), l.split()[4]) for l in open(f'{O}/regression819.txt') if l.startswith('old vs'))
prb = dict(((l.split()[1], l.split()[4]), l.split()[5]) for l in open(f'{O}/probes819.txt') if l.startswith('probe'))
okO = all(reg[('new_off', t)] == 'SAME' for t in ('1', '2', 'L2', '5Ps', '5H', '6')) and all(prb[(t, 'new_off')] == 'SAME' for t in ('H2p', 'H3p', 'H4p', 'H3s', 'H3d', 'Ps2p'))
pr(f"H819o (new_off = old: regresja 6/6, sondy 6/6): {'TAK' if okO else 'NIE'}  " + ' '.join(f'{t}:{reg[("new_off", t)]}' for t in ('1', '2', 'L2', '5Ps', '5H', '6')))
okH = reg[('new', '5H')] == 'SAME' and all(prb[(t, 'new')] == 'SAME' for t in ('H2p', 'H3p', 'H4p', 'H3s', 'H3d'))
pr(f"H819h (wodor new = old: 5H i sondy H2p, H3p, H4p, H3s, H3d): {'TAK' if okH else 'NIE'}  " + ' '.join(f'{t}:{prb[(t, "new")]}' for t in ('H2p', 'H3p', 'H4p', 'H3s', 'H3d', 'Ps2p')))
pr('  RAPORT zestawy Ps new vs old: ' + ' '.join(f'{t}:{reg[("new", t)]}' for t in ('1', '2', 'L2', '5Ps', '6')))
# ---------------- H819v
V = {v: open(f'{O}/validation_{v}.txt').read() for v in ('old', 'new')}
okV = 'EXIT 0' in V['new'] and 'Enforced checks:     72/72' in V['new']
pr(f"H819v (walidacja new EXIT 0, 72/72): {'TAK' if okV else 'NIE'}  (old: {'EXIT 0' in V['old']} {'72/72' in V['old']})")
dl = [(a, b) for a, b in zip(V['old'].splitlines(), V['new'].splitlines()) if a != b]
pr(f'  RAPORT walidacja: {len(dl)} roznych wierszy old / new' + ''.join(f'\n    old: {a[:150]}\n    new: {b[:150]}' for a, b in dl[:12]))
# ---------------- H819s (production statistics; parser stat819, gate on the old files)
import stat819 as S9
miss = S9.check({t: f'{O}/reg/old_{t}.out' for t in S9.PAT}); VALID = not miss
pr(f"BRAMKA PARSERA (pliki old 819): {'OK' if VALID else 'BRAK ' + str(miss) + ' -> AUDYT NIEWAZNY'}")
okS = VALID
for t in S9.PAT:
    so, sn = S9.parse(f'{O}/reg/old_{t}.out', t), S9.parse(f'{O}/reg/new_{t}.out', t)
    for k, _, kind in S9.PAT[t]:
        if k not in so or k not in sn: okS = False; pr(f'  {t} {k}: BRAK wiersza'); continue
        vo, vn = so[k][1], sn[k][1]
        if kind == 'pm':
            (a_, ea), (b_, eb) = vo, vn; ok = abs(b_ - a_) <= 2*ea
            pr(f'  zestaw {t} {k}: old {a_:.6g} +/- {ea:.4g}, new {b_:.6g} +/- {eb:.4g}; (new - old)/sigma_old = {(b_ - a_)/ea:+.3f} {"" if ok else "<-- poza 2 sigma"}')
        elif kind == 'eq':
            ok = vo == vn; pr(f'  zestaw {t} {k}: old {vo[0]:g}, new {vn[0]:g} {"" if ok else "<-- rozne"}')
        else:
            ok = vo[0] == 0 and vn[0] == 0; pr(f'  zestaw {t} {k}: old {vo[0]:g}, new {vn[0]:g} {"" if ok else "<-- nie zero"}')
        okS &= ok
pr(f"H819s (GLOWNA 2: statystyka Ps new vs old w 2 sigma_old, ukonczenie rowne, 0 failed): {'TAK' if okS else 'NIE'}")
# ---------------- H819c (wall time)
W = {}
for l in open(f'{O}/walltime819.txt'):
    v, t, s = l.split(); W[(v, t)] = float(s)
okC = True
for t in ('1', '2', 'L2', '5Ps', '6'):
    r = W[('new', t)]/max(W[('old', t)], 1.0); okC &= r <= 1.5
    pr(f'  czas {t}: old {W[("old", t)]:.0f} s, new {W[("new", t)]:.0f} s, new_off {W[("new_off", t)]:.0f} s; new/old {r:.3f}')
pr(f"H819c (koszt: new/old <= 1,5 w kazdym zestawie Ps): {'TAK' if okC else 'NIE'}")
# ---------------- H819r (orbit n = 2 under the new default)
x, w = np.polynomial.legendre.leggauss(5); u = 0.75 + 0.25*x[::-1]
NODES = dict(zip([f'u5{c}' for c in 'abcde'], zip(1/u, 0.25*w[::-1])))
nu_of = lambda Om: sum(w*L**2*Om[k] for k, (L, w) in NODES.items())/(2*math.pi)
def read(path, tag='mechB0'):   # = ana815.read
    rows = []; has_sum = has_cons = False
    for ln in open(path):
        if ln.startswith(f'PAIR TORQUE799 {tag} sum'): has_sum = True
        if ln.startswith(f'PAIR CONS814 {tag} start'): has_cons = True
        if not ln.startswith(f'PAIR APSE808 {tag} '): continue
        w_ = ln.split(); tc = float(w_[w_.index('tc') + 1])
        r = [np.array([float(v) for v in w_[w_.index(f'r{q}') + 1:w_.index(f'r{q}') + 4]]) for q in range(3)]
        rows.append((tc, r))
    r0 = rows[0][1]; Lh = np.cross(r0[0], r0[2]); Lh /= np.linalg.norm(Lh)
    e1 = r0[1]/np.linalg.norm(r0[1]); e1 -= Lh*np.dot(Lh, e1); e1 /= np.linalg.norm(e1); e2 = np.cross(Lh, e1)
    T = []; Pp = []; nbad = 0; prev = None
    for tc, r in rows:
        th = np.unwrap([math.atan2(np.dot(v, e2), np.dot(v, e1)) for v in r]); rho = [np.linalg.norm(v) for v in r]
        xx = th - th[1]; Am = np.vstack([np.ones(3), xx, xx**2]).T; a0, b, cq = np.linalg.solve(Am, rho)
        xs = -b/(2*cq); hh = max(abs(xx[0]), abs(xx[2]))
        if not (cq < 0 and xx[0] - 0.5*hh <= xs <= xx[2] + 0.5*hh): nbad += 1
        phi = th[1] + xs
        if prev is not None: phi += 2*math.pi*round((prev - phi)/(2*math.pi))
        prev = phi; T.append(tc); Pp.append(phi)
    return dict(lin=np.polyfit(np.array(T), np.array(Pp), 1)[0], n=len(T), nbad=nbad, has_sum=has_sum, has_cons=has_cons, path=path)
A = {(s, k): read(f'{O}/runs/N/F1_ops_{k}_s{s}_5.err') for s in (2048, 4096) for k in NODES}
PP = {k: read(f'{O}/runs/PN/F0_pps_{k}_s2048_5.err') for k in NODES}
allf = list(A.values()) + list(PP.values())
okF = len(allf) == 15 and all(v['has_sum'] and v['has_cons'] and v['n'] >= 195 and v['nbad'] == 0 for v in allf)
pr(f'  pliki orbity: {len(allf)}, min apocentrow {min(v["n"] for v in allf)}, zlych {sum(v["nbad"] for v in allf)} -> {"OK" if okF else "BRAK"}')
okR = okF; R2, Pe, P0 = {}, {}, {}
for k in NODES:
    y = [A[(s, k)]['lin'] for s in (2048, 4096)]; R2[k] = y[1] + (y[1] - y[0])/3
    st = parse_start(A[(4096, k)]['path']); Pe[k] = P15.omega(st, 0.01)[0]; P0[k] = P15.omega(st, 0.0)[0]
    rr = R2[k]/Pe[k] - 1; okR &= abs(rr) <= 5e-6
    pr(f'  {k}: Omega 2048/4096 {y[0]:.10e} / {y[1]:.10e}, R2 {R2[k]:.10e}; pred(0,01) {Pe[k]:.10e} ({rr:+.2e}); pred(0) {P0[k]:.10e}')
n2, ne, n0 = nu_of(R2), nu_of(Pe), nu_of(P0)
dc = max(abs(PP[k]['lin']/A[(2048, k)]['lin'] - 1) for k in NODES)
okR &= abs(n2 - ne) <= 0.010e6 and abs(n2 - n0) <= 0.070e6 and dc <= 1e-8
pr(f'  o-Ps (PHEN=2) nu_orb R2 = {units(n2)}; nu_pred(0,01) {ne/1e6:.4f} ({(n2 - ne)/1e6:+.4f}); nu_pred(0) {n0/1e6:.4f}; oslona w produkcji {(n2 - n0)/1e6:+.4f} MHz')
nP = nu_of({k: PP[k]['lin'] for k in NODES}) + (n2 - nu_of({k: A[(2048, k)]['lin'] for k in NODES}))
pr(f'  p-Ps (PHEN=1) 2048 + poprawka kroku z o-Ps: nu_orb = {units(nP)}; p-Ps / o-Ps max {dc:.2e}')
pr(f"H819r (GLOWNA 1: orbita = pred(0,01) do 0,010 MHz i 5e-6, oslona <= 0,070 MHz, p-Ps = o-Ps 1e-8): {'TAK' if okR else 'NIE'}")
import tarfile, tempfile, filecmp
tmp = tempfile.mkdtemp(); tarfile.open(os.path.join(D, '..', 'section-818', 'runs', 'runs818.tar.xz')).extractall(tmp)
same = [filecmp.cmp(v['path'], os.path.join(tmp, os.path.relpath(v['path'], O)), shallow=False) for v in allf]
pr(f'  RAPORT orbita 819 vs 818 bajt w bajt: {sum(same)}/{len(same)}')
allok = okO and okH and okV and okS and okC and okR
pr(f"DECYZJA: {'AUDYT NIEWAZNY (parser)' if not VALID else 'domyslna skala oslony 0,01 r* PRZYJETA' if allok else 'zmiana COFNIETA, domyslna 0,05 zostaje'}")

# Audit 699: start increments of L^, s1^, s2^ over N orbits from the 2P formation state -- secular transport, mechanical path (and its precession
# sectors), and the mean field (Hartree) of the product state under H_2P and its parts (model QM reference, not a measurement).
# usage: python3 ana699.py <699 runs dir with def/ charge/ dipole/ ctl/> <698 out/post698.txt>
import sys, glob, math, os, re
import numpy as np
HBAR = 1.054571817e-34
PAIRS = (('antiproton', 'antihydrogen'), ('proton', 'hydrogen'))
FIELDS = ('1.0329', '0.1', '0.01')
R2 = math.sqrt(0.5)
LX = np.array([[0, R2, 0], [R2, 0, R2], [0, R2, 0]], complex); LD = np.array([[0, R2, 0], [-R2, 0, R2], [0, -R2, 0]], complex)
LZ = np.diag([1.0, 0.0, -1.0]).astype(complex)
SX = np.array([[0, .5], [.5, 0]], complex); SD = np.array([[0, .5], [-.5, 0]], complex); SZ = np.diag([.5, -.5]).astype(complex)
LOPS = (LX, -1j*LD, LZ); SOPS = (2*SX, -2j*SD, 2*SZ)   # <L> and <sigma>


def cohL(u):   # psCoherentL1
    ct = max(-1.0, min(1.0, u[2])); ph = math.atan2(u[1], u[0]); e = complex(math.cos(ph), -math.sin(ph))
    return np.array([e*(1 + ct)*0.5, math.sqrt(max(0.0, 1 - ct*ct))/math.sqrt(2.0), e.conjugate()*(1 - ct)*0.5])


def cohS(u):   # psCoherentHalf
    th = math.acos(max(-1.0, min(1.0, u[2]))); ph = math.atan2(u[1], u[0])
    return np.array([np.exp(-0.5j*ph)*math.cos(0.5*th), np.exp(0.5j*ph)*math.sin(0.5*th)])


def prod(l, a, b): return np.einsum('i,p,q->ipq', l, a, b).reshape(12)


def ev(v, ops): return np.array([np.real(np.vdot(v, o @ v)) for o in ops])


def parse(path):
    d = {'H': {}, 'inc': {}}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR QM1S699 H'): d['H'][t[3]] = np.array([float(x) for x in t[4:148]]).reshape(12, 12)
        elif line.startswith('PAIR QM1S699 '):
            g = {t[i]: np.array([float(x) for x in t[i+1:i+4]]) for i in range(len(t)) if t[i] in ('L0', 's10', 's20', 'L1', 's11', 's21')}
            g['dt'] = float(t[t.index('dt')+1]); g['Ln0'] = float(t[t.index('Ln0')+1]); d['inc'][t[2]] = g
        elif line.startswith('PAIR QM1S695 form'):
            d['w'] = np.array([float(x) for x in t[4:16]]); d['p'] = np.array([float(x) for x in t[17:29]])
    return d


def hartree(H, X0, dt, steps=400):
    l, a, b = cohL(X0[0]), cohS(X0[1]), cohS(X0[2]); T = H.reshape(3, 2, 2, 3, 2, 2)
    def rhs(l, a, b):
        hl = np.einsum('p,q,ipqjrs,r,s->ij', a.conj(), b.conj(), T, a, b)
        ha = np.einsum('i,q,ipqjrs,j,s->pr', l.conj(), b.conj(), T, l, b)
        hb = np.einsum('i,p,ipqjrs,j,r->qs', l.conj(), a.conj(), T, l, a)
        return -1j/HBAR*(hl @ l), -1j/HBAR*(ha @ a), -1j/HBAR*(hb @ b)
    h = dt/steps
    for _ in range(steps):
        k1 = rhs(l, a, b); k2 = rhs(*(x + 0.5*h*k for x, k in zip((l, a, b), k1)))
        k3 = rhs(*(x + 0.5*h*k for x, k in zip((l, a, b), k2))); k4 = rhs(*(x + h*k for x, k in zip((l, a, b), k3)))
        l, a, b = (x + h/6*(p + 2*q + 2*r + s) for x, p, q, r, s in zip((l, a, b), k1, k2, k3, k4))
        l, a, b = l/np.linalg.norm(l), a/np.linalg.norm(a), b/np.linalg.norm(b)
    return [ev(l, LOPS), ev(a, SOPS), ev(b, SOPS)]


def unit(v): return v/np.linalg.norm(v)


def pops(V, X): return np.abs(V.T @ prod(cohL(unit(X[0])), cohS(unit(X[1])), cohS(unit(X[2]))))**2


def main():
    R, P8 = sys.argv[1:3]
    c698 = {(m.group(1), m.group(2), m.group(3)): float(m.group(4)) for m in
            re.finditer(r'(\w+)\s+B\s+([\d.]+) T seed (\d): c ([\d.]+)', open(P8).read())}
    out = []
    # operator convention check: <L>, <sigma> of the coherent states = their axis
    u = unit(np.array([0.3, -0.5, 0.7])); chk = max(np.abs(ev(cohL(u), LOPS) - u).max(), np.abs(ev(cohS(u), SOPS) - u).max())
    out.append(f'convention: |<L> - u|, |<sigma> - u| for a coherent state {chk:.1e}')
    fs = sorted(glob.glob(f'{R}/def/*.err'))
    same = sum(''.join(l for l in open(f) if not l.startswith('PAIR QM1S699')) == open(f'{R}/ctl/{os.path.basename(f)}').read() for f in fs)
    out.append(f'H699a (files): default sector with CREM_PAIR_TORQUE without QM1S699 = without it {same}/{len(fs)}: {same == len(fs) == 12}')
    okB = okC = True; okD = None; okE = True; culprit = {}
    for tag, name in PAIRS:
        for B in FIELDS:
            for s in ('0', '1'):
                key = f'{tag}_B{B}_A2_{s}.err'; d = parse(f'{R}/def/{key}'); H = d['H']
                w, V = np.linalg.eigh(H['all']); inc = d['inc']; g = inc['sec']
                X0 = [g['L0'], g['s10'], g['s20']]; dt = g['dt']
                p0 = pops(V, X0); pf = d['p'][np.argsort(d['w'])]
                hsum = np.abs(H['fs'] + H['hfs'] + H['z'] - H['all']).max()/np.abs(H['all']).max()
                okB &= np.abs(p0 - pf).max() <= 1e-6 and hsum <= 1e-9
                paths = {k: [inc[k]['L1'], inc[k]['s11'], inc[k]['s21']] for k in ('sec', 'secB0', 'mech', 'mechB0')}
                for sec in ('charge', 'dipole'):
                    e = parse(f'{R}/{sec}/{key}')['inc']
                    if 'mechB0' in e: paths[f'{sec}B0'] = [e['mechB0']['L1'], e['mechB0']['s11'], e['mechB0']['s21']]
                for part, Hp in (('Hall', H['all']), ('HB0', H['fs'] + H['hfs']), ('Hfs', H['fs']), ('Hhfs', H['hfs'])):
                    paths[part] = hartree(Hp, X0, dt)
                dX = {k: [np.asarray(v[i]) - X0[i] for i in range(3)] for k, v in paths.items()}
                dD = {k: pops(V, v) - p0 for k, v in paths.items()}
                c = float(dD['mech'] @ dD['sec']/(dD['sec'] @ dD['sec'])); cr = c698.get((name, B, s), float('nan'))
                okC &= abs(c/cr - 1) <= 0.25
                rel = lambda a, b: float(np.linalg.norm(a - b)/np.linalg.norm(b)) if np.linalg.norm(b) > 0 else float('nan')
                zs = [dX['sec'][i] - dX['secB0'][i] for i in range(3)]; zm = [dX['mech'][i] - dX['mechB0'][i] for i in range(3)]
                zr = [rel(zs[i], zm[i]) for i in range(3)]; ir = [rel(dX['secB0'][i], dX['mechB0'][i]) for i in range(3)]
                if name == 'hydrogen' and B == '1.0329': okD = (okD is not False) and max(zr) > 0.10
                if name == 'hydrogen' and B == '0.01': okE &= max(ir) > 0.10
                # contribution of each vector's (sec - mech) difference to |dD_sec - dD_mech|
                contrib = []
                for i in range(3):
                    Xs = [paths['mech'][j] if j != i else paths['sec'][j] for j in range(3)]
                    contrib.append(float(np.linalg.norm(pops(V, Xs) - pops(V, paths['mech']))))
                culprit[(name, B, s)] = ('L', 's1', 's2')[int(np.argmax(contrib))]
                out.append(f'{name:12s} B {B:>6} T seed {s} dt {dt*1e12:.4f} ps | start p = form {np.abs(p0 - pf).max():.1e}, |fs+hfs+z-all| {hsum:.1e}'
                           f' | c699 {c:.3f} (698 {cr:.3f}) | Ln0 {g["Ln0"]:.4f}')
                for i, v in enumerate(('L', 's1', 's2')):
                    cells = ' '.join(f'{k} {np.linalg.norm(dX[k][i]):.3e}' for k in dX)
                    out.append(f'   |d{v}^| {cells}')
                    out.append(f'      rel|sec - mech| {rel(dX["sec"][i], dX["mech"][i]):.3f}  Zeeman part rel {zr[i]:.3f} (|sec| {np.linalg.norm(zs[i]):.2e} |mech| {np.linalg.norm(zm[i]):.2e})'
                               f'  B0 part rel {ir[i]:.3f}  | vs Hartree all: sec {rel(dX["sec"][i], dX["Hall"][i]):.3f} mech {rel(dX["mech"][i], dX["Hall"][i]):.3f}'
                               + (f'  | charge B0 vs Hfs {rel(dX["chargeB0"][i], dX["Hfs"][i]):.3f}  dipole B0 vs Hhfs {rel(dX["dipoleB0"][i], dX["Hhfs"][i]):.3f}'
                                  if 'chargeB0' in dX and 'dipoleB0' in dX else ''))
                out.append('   |dD| ' + ' '.join(f'{k} {np.linalg.norm(dD[k]):.2e}' for k in dD) + ' | contribution of (sec - mech) per vector L s1 s2 '
                           + ' '.join(f'{x:.2e}' for x in contrib) + f' -> {culprit[(name, B, s)]}')
    okF = all(culprit[(n, B, '0')] == culprit[(n, B, '1')] for _, n in PAIRS for B in FIELDS)
    out.append(f'H699b (start populations = form to 1e-6, parts sum to all to 1e-9): {okB}')
    out.append(f'H699c (c699 within 25 % of c698, all 12): {okC}')
    out.append(f'H699d (hydrogen 1 T: Zeeman part sec vs mech differs > 10 % for some vector, both seeds): {okD}')
    out.append(f'H699e (hydrogen 0.01 T: B0 part sec vs mech differs > 10 % for some vector, both seeds): {okE}')
    out.append(f'H699f (same culprit vector in both seeds of every cell): {okF}  ' + ', '.join(f'{n} {B}: {culprit[(n, B, "0")]}/{culprit[(n, B, "1")]}' for _, n in PAIRS for B in FIELDS))
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana699.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

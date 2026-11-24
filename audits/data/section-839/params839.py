# Audit 839: orbit (r_p, e from START838) and dipoles (CREM_LEDGER_START) in the frame P = periapsis, Q, Z = L, from the 838 ZA files.
import sys, numpy as np
R = sys.argv[1]
SYS = {'ops': ('F1_ops', 'electron,positron'), 'H0': ('F0_H', 'proton,electron'), 'H1': ('F1_H', 'proton,electron')}
for sk, (nm, pair) in SYS.items():
    for k in 'abcde':
        L = open(f'{R}/ZA/D_{nm}_u5{k}_s2048.err').read().splitlines()
        led = [l for l in L if l.startswith('CREM_LEDGER_START')][0].split()
        m1 = np.array(list(map(float, led[2:5]))); m2 = np.array(list(map(float, led[6:9])))
        c = [l for l in L if l.startswith('PAIR CONS814')][0].split(); i = c.index('period'); v = list(map(float, c[i+2:i+14]))
        x1, x2, v1, v2 = [np.array(v[3*q:3*q+3]) for q in range(4)]
        r = x1 - x2; vv = v1 - v2; Z = np.cross(r, vv); Z /= np.linalg.norm(Z); P = -r/np.linalg.norm(r); Q = np.cross(Z, P)
        st = [l for l in L if 'START838' in l][0].split()
        f = lambda m: [m @ P, m @ Q, m @ Z]
        print(sk, f'u5{k}', pair, st[st.index('rp')+1], st[st.index('e')+1], *[f'{x:.15e}' for x in f(m1) + f(m2)])

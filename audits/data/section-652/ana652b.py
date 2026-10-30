# Audit 652, post hoc (recorded in supplement 2): per-seed split of N = jump - classical.
# N = (Qa - Qb) + (Zb - Za) + R',  Qb = Q_form + T (transport: dU + dZ in 2P);  so
# N = -Q_form - T + Qa + (Zb - Za) + R'.  Here each term per channel, and the rms of N left after each step.
import sys,math; sys.argv=[sys.argv[0],sys.argv[1]]
exec(open(__file__.replace('ana652b.py','ana652.py')).read().split('def stat')[0])
rms=lambda x: math.sqrt(sum(v*v for v in x)/len(x))
for B in ('1e-6','0.1','1'):
    for c,name in (('O1','o-Ps'),('P1','p-Ps')):
        R=[parse(f'{d}/B{B}_{c}_L2_{i}.err') for i in range(16)]
        rows=[]
        for r in R:
            p=r['ph']; e=r['eig']; N=(e['jump']-e['cl'])*u; Qf=r['form'][3]*u
            T=p['Qb']*u-Qf; Qa=p['Qa']*u; dZ=(p['Zb']-p['Za'])*u
            Rp=N-(-Qf-T+Qa+dZ); rows.append((N,Qf,T,Qa,dZ,Rp,p['dorb']*u))
        N=[x[0] for x in rows]
        r1=[x[0]+x[1] for x in rows]; r2=[x[0]+x[1]+x[2] for x in rows]; r3=[x[0]+x[1]+x[2]-x[3] for x in rows]; r4=[x[0]+x[1]+x[2]-x[3]-x[4] for x in rows]
        print(f'{name} B {B} T: rms N {rms(N):7.3f} | after -Q_form {rms(r1):7.3f} | after transport {rms(r2):7.3f} | after Q_after {rms(r3):7.3f} | after Zeeman {rms(r4):.2e}  (R\' max {max(abs(x[5]) for x in rows):.1e}, rung delta rms {rms([x[6] for x in rows]):.3f})')
        print(f'      rms terms: Q_form {rms([x[1] for x in rows]):7.3f}  transport T {rms([x[2] for x in rows]):7.3f}  Q_after {rms([x[3] for x in rows]):7.3f}  Zb-Za {rms([x[4] for x in rows]):7.3f}  ueV')

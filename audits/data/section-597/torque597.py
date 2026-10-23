# Audit 597: the classical radiation-reaction torque never reverses the orbit; E1 channels from "no reversal" + ellipse existence.
import numpy as np
for e in (0.1,0.66,0.95,0.999):
    a=1.0; M=np.linspace(0,2*np.pi,4096,endpoint=False); u=np.full_like(M,np.pi)
    for _ in range(200): u-=(u-e*np.sin(u)-M)/(1-e*np.cos(u))
    r=np.stack([a*(np.cos(u)-e), a*np.sqrt(1-e*e)*np.sin(u), 0*u],1)
    n=1.0; udot=n/(1-e*np.cos(u))                                   # dE/dt, mean motion 1 (k = mu = a = 1)
    v=np.stack([-a*np.sin(u)*udot, a*np.sqrt(1-e*e)*np.cos(u)*udot, 0*u],1)
    R=np.linalg.norm(r,axis=1)[:,None]
    adot=-v/R**3+3*r*np.sum(r*v,1)[:,None]/R**5                   # d/dt(-r/r^3)
    L=np.cross(r,v); tq=np.cross(r,adot); ref=-L/R**3
    dev=np.max(np.linalg.norm(tq-ref,axis=1)/np.linalg.norm(ref,axis=1))
    print(f"T1  e={e}: max |r x a' + L/r^3| / |L/r^3| = {dev:.1e};  L_z sign constant: {bool(np.all(L[:,2]>0))}")
mis=0; extraA=0; tot=0
for n in range(2,31):
    for l in range(n):
        L=l+0.5
        for n2 in range(1,n):
            for l2 in (l-1,l+1):
                if l2<0: continue
                tot+=1
                co=(l2==l-1)
                crit=(L>=1.0 if co else True) and (l2+0.5)<=n2
                qm=l2<=n2-1
                mis+=crit!=qm
            # vector rule A': co-rotating from L = 1/2 gives |L'| = 1/2 (l' = l = 0): s -> s
            if l==0 and n2>=1: extraA+=1
print(f"T2  channels {tot}: criterion 'co-rotating only at |L| >= hbar' + ellipse vs QM: {mis} mismatches;  vector rule A' adds {extraA} s -> s channels (Delta l = 0)")
# T3: lattice offsets before/after each engine photon (modes592.txt)
offs=[]
for line in open('../section-592/modes592.txt'):
    t=line.split()
    if len(t)<6 or not t[0].endswith('.err'): continue
    Lb,La=float(t[2]),float(t[3]); f=lambda x: abs((x-0.5)-round(x-0.5))
    offs.append((f(Lb),f(La)))
print(f"T3  {len(offs)} photons: max |L/hbar - 1/2 - nearest integer| before {max(o[0] for o in offs):.1e}, after {max(o[1] for o in offs):.1e} (printed L to 5-6 digits)")

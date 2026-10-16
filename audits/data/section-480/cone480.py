# Audit 480: spin on a cone with a fast phase vs the hbar/2 vector -- averages of all terms linear in S, and the non-linear ones.
import math, random
random.seed(480)
hbar=1.0; c=1.0
def add(a,b): return [a[i]+b[i] for i in range(3)]
def sub(a,b): return [a[i]-b[i] for i in range(3)]
def mul(a,k): return [x*k for x in a]
def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def cross(a,b): return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def norm(a): return math.sqrt(dot(a,a))
def unit(a): return mul(a,1/norm(a))
def rnd(): return [random.gauss(0,1) for _ in range(3)]
def basis(a):
    e1=unit(cross(a,[1,0,0] if abs(a[0])<0.9 else [0,1,0])); e2=cross(a,e1); return e1,e2
t=hbar/math.sqrt(2); proj=hbar/2
def cone(a,psi):
    e1,e2=basis(a); return add(mul(a,proj),add(mul(e1,t*math.cos(psi)),mul(e2,t*math.sin(psi))))
def dipB(R,m):   # field of a point magnetic dipole m at displacement R (units mu0/4pi = 1)
    r=norm(R); n=mul(R,1/r); return mul(sub(mul(n,3*dot(n,m)),m),1/r**3)
N=4096; worst={}
def rec(k,v): worst[k]=max(worst.get(k,0.0),v)
for trial in range(200):
    a=unit(rnd()); gam=random.uniform(-3,3); Om=rnd(); B=rnd(); E=rnd(); L=mul(rnd(),2.0); A=random.uniform(0.1,2); R=rnd(); v=mul(rnd(),0.1)
    S0=mul(a,proj)
    av={k:[0,0,0] for k in ('torque','J','react','Phid','Jhid','Bdip','dJ')}; avs={'E':0.0,'SO':0.0}
    for i in range(N):
        psi=2*math.pi*i/N; S=cone(a,psi); mu=mul(S,gam)
        Beff=sub(B,mul(cross(v,E),0.5))
        terms={'torque':cross(Om,S),'J':add(L,S),'react':mul(cross(Om,S),-1),'Phid':cross(mu,E),'Jhid':cross(R,cross(mu,E)),'Bdip':dipB(R,mu),
               'dJ':add(mul(cross(Om,S),-1),cross(Om,S))}
        for k,vv in terms.items(): av[k]=add(av[k],mul(vv,1/N))
        avs['E']+=-dot(mu,Beff)/N; avs['SO']+=A*dot(L,S)/N
    mu0=mul(S0,gam); Beff=sub(B,mul(cross(v,E),0.5))
    ref={'torque':cross(Om,S0),'J':add(L,S0),'react':mul(cross(Om,S0),-1),'Phid':cross(mu0,E),'Jhid':cross(R,cross(mu0,E)),'Bdip':dipB(R,mu0),'dJ':[0,0,0]}
    for k in av:
        sc=max(norm(ref[k]),norm(terms[k]),1e-300); rec('(lin) '+k,norm(sub(av[k],ref[k]))/sc)
    rec('(lin) energy -mu.B',abs(avs['E']+dot(mu0,Beff))/max(abs(dot(mu0,Beff)),1e-300)); rec('(lin) spin-orbit A L.S',abs(avs['SO']-A*dot(L,S0))/max(abs(A*dot(L,S0)),1e-300))
print('H480a: worst relative difference (cone average vs hbar/2 vector), 200 random configurations, 4096 phase nodes:')
for k in sorted(worst): print(f'   {k:28s} {worst[k]:.2e}')
# (10) S1.S2 of two cones on a common axis, opposite projections (T0 / singlet) and same projections
print('\nH480b: <S1.S2> on a common axis (hbar^2):')
a=[0,0,1]
for name,p1,p2 in (('opposite projections',+1,-1),('same projections',+1,+1)):
    for ph,label in ((None,'free phases'),(0.0,'locked, in phase'),(math.pi,'locked, antiphase')):
        s=0.0; M=512
        for i in range(M):
            for j in range(M if ph is None else 1):
                psi1=2*math.pi*i/M; psi2=2*math.pi*j/M if ph is None else psi1+ph
                S1=cone(a,psi1); S1[2]*=p1; S2=cone(a,psi2); S2[2]*=p2
                s+=dot(S1,S2)/(M*(M if ph is None else 1))
        print(f'   {name:22s} {label:18s}: {s:+.6f}   (hbar/2 vectors: {p1*p2*0.25:+.3f})')
# (11) |S|^2, |J|^2
L=[0,0,1.5]; aS=unit([1,1,0.3]); s2=0; j2=0; M=4096
for i in range(M):
    S=cone(aS,2*math.pi*i/M); s2+=dot(S,S)/M; j2+=dot(add(L,S),add(L,S))/M
S0=mul(aS,proj)
print(f'\nH480c: <|S|^2> = {s2:.6f} (hbar/2: {dot(S0,S0):.6f}); <|J|^2> = {j2:.6f} (hbar/2: {dot(add(L,S0),add(L,S0)):.6f}); difference {j2-dot(add(L,S0),add(L,S0)):.6f} = t^2 = {t*t:.6f}')
# (12) M1 power ~ |mu''|^2 with fast phase omega_f vs slow precession Omega (hbar/2: |Omega x (Omega x mu)|^2)
wf=1e3; Om=[0,0,1.0]; gam=1.0; a=unit([0.3,0.2,1]); acc=0; M=4096
for i in range(M):
    psi=2*math.pi*i/M; e1,e2=basis(a); mdd=mul(add(mul(e1,math.cos(psi)),mul(e2,math.sin(psi))),-t*gam*wf*wf); acc+=dot(mdd,mdd)/M
m0=mul(a,proj*gam); ref=cross(Om,cross(Om,m0))
print(f'H480d: <|mu\'\'|^2> fast phase (omega_f = {wf:g} Omega) {acc:.4e}  vs hbar/2 vector {dot(ref,ref):.4e}  (ratio {acc/dot(ref,ref):.3e})')
# (13) instantaneous J with a free fast phase: dS/dt = omega_f a x S has no source torque
print(f'       free fast phase: |dS/dt| = omega_f t = {wf*t:.4e} (in units of hbar Omega) with no torque -- J = L + S oscillates with amplitude t = {t:.4f} hbar; its average over the phase stays L + S0')

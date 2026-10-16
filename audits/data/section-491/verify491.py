# Audit 491: verification of the 3P_J claims (lock of the spin transverse phase on the position azimuth).
import math, random
random.seed(491)
# (1) T1: <n n r^-3>_t over a Kepler orbit, in the orbital plane
print('(1) <n_x^2 r^-3>/<r^-3>, <n_y^2 r^-3>/<r^-3>, <n_x n_y r^-3>/<r^-3> (periapsis along x), time average:')
for e in (0.0,0.5,0.9,0.99):
    N=400000; sxx=syy=sxy=w=0
    for i in range(N):
        E=2*math.pi*(i+0.5)/N; dt=1-e*math.cos(E); r=1-e*math.cos(E)
        x=(math.cos(E)-e)/r; y=math.sqrt(1-e*e)*math.sin(E)/r; f=dt/r**3
        sxx+=f*x*x; syy+=f*y*y; sxy+=f*x*y; w+=f
    print(f'    e {e:.2f}: {sxx/w:.6f} {syy/w:.6f} {sxy/w:+.1e}')
# (2) model with L as a cone
def cg(j1,m1,j2,m2,J,M):
    if m1+m2!=M or abs(m1)>j1 or abs(m2)>j2 or abs(M)>J: return 0.0
    f=math.factorial
    pre=math.sqrt((2*J+1)*f(J+j1-j2)*f(J-j1+j2)*f(j1+j2-J)/f(j1+j2+J+1)*f(J+M)*f(J-M)*f(j1-m1)*f(j1+m1)*f(j2-m2)*f(j2+m2))
    s=0
    for k in range(20):
        a=[k,j1+j2-J-k,j1-m1-k,j2+m2-k,J-j2+m1+k,J-j1-m2+k]
        if min(a)<0: continue
        s+=(-1)**k/(f(a[0])*f(a[1])*f(a[2])*f(a[3])*f(a[4])*f(a[5]))
    return pre*s
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def unit(a): n=math.sqrt(dot(a,a)); return tuple(x/n for x in a)
def S12_J(J,Llen,slen,c,delta,lockL=None,NS=6000):
    M=0 if J==0 else J
    comps=[(ml,ms,cg(1,ml,1,ms,J,M)**2) for ml in (-1,0,1) for ms in (-1,0,1) if ml+ms==M]
    t=math.sqrt(slen*slen-1); tot=0
    for ml,ms,w in comps:
        if w<1e-12: continue
        ct=ml/Llen; st=math.sqrt(max(0,1-ct*ct)); acc=0
        for k in range(NS):
            chi=2*math.pi*random.random() if lockL is None else lockL
            Lh=(st*math.cos(chi),st*math.sin(chi),ct)
            ref=(0,0,1) if abs(ct)<0.99 else (1,0,0)
            e1=unit(cross(ref,Lh)); e2=cross(Lh,e1); ph=2*math.pi*random.random()
            n=tuple(math.cos(ph)*e1[i]+math.sin(ph)*e2[i] for i in range(3))
            al=math.atan2(n[1],n[0])
            p1,p2=(1,1) if ms==1 else ((-1,-1) if ms==-1 else (1,-1))
            if c is None: ps1=2*math.pi*random.random(); ps2=ps1+delta if (ms==0 and delta is not None) else 2*math.pi*random.random()
            else: ps1=al+c; ps2=ps1+delta if (delta is not None) else 2*math.pi*random.random()
            s1=(t*math.cos(ps1),t*math.sin(ps1),p1); s2=(t*math.cos(ps2),t*math.sin(ps2),p2)
            acc+=3*dot(s1,n)*dot(s2,n)-dot(s1,s2)
        tot+=w*acc/NS
    return tot
for label,Llen,slen in (('QM lengths |L| = sqrt2, |sigma| = sqrt3',math.sqrt(2),math.sqrt(3)),('Langer |L| = 1.5, |sigma| = 2',1.5,2.0)):
    print(f'\n(2) {label}')
    for J in (0,1,2):
        free=S12_J(J,Llen,slen,None,0.0); free_pi=S12_J(J,Llen,slen,None,math.pi)
        best_min=(1e9,None); best_max=(-1e9,None)
        for dl in (0.0,math.pi):
            for k in range(24):
                c=2*math.pi*k/24; v=S12_J(J,Llen,slen,c,dl,NS=2500)
                if v<best_min[0]: best_min=(v,(math.degrees(c),dl))
                if v>best_max[0]: best_max=(v,(math.degrees(c),dl))
        print(f'   J={J}: no lock (T0 in phase) {free:+.3f}, no lock (T0 antiphase) {free_pi:+.3f}; lock scan min {best_min[0]:+.3f} at c={best_min[1][0]:.0f} deg delta={best_min[1][1]:.2f}, max {best_max[0]:+.3f} at c={best_max[1][0]:.0f} deg delta={best_max[1][1]:.2f}   (QM {[-4,2,-0.4][J]:+.2f})')

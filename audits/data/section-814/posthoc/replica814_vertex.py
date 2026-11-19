# POST HOC 814: replica of the engine's mechanical step (crem_trajectory/electrodynamics) with conservative forces.
import sys, math, numpy as np
sys.path.insert(0,'.')
from pred814 import parse_start, darwin_on_first
f=sys.argv[1]; h=float(sys.argv[2]); orbits=int(sys.argv[3]); cap=float(sys.argv[4]) if len(sys.argv)>4 else 5e-18
st=parse_start(f); m1,m2,kq,c=st['m1'],st['m2'],st['kq'],st['c']; K=-kq; mu=m1*m2/(m1+m2)
def vel(p,m): return p/(m*math.sqrt(1+(p@p)/(m*m*c*c)))
def mom(v,m): return v*m/math.sqrt(1-(v@v)/(c*c))
def forces(x1,x2,v1,v2):
    d=x1-x2; r=math.sqrt(d@d); F1=d*(kq/r**3); F2=-F1
    return F1+darwin_on_first(v1,v2,F2/m2,d,kq,c), F2+darwin_on_first(v2,v1,F1/m1,-d,kq,c)
x1,x2,v1,v2=st['x1'].copy(),st['x2'].copy(),st['v1'].copy(),st['v2'].copy()
d0=x1-x2; Lh=np.cross(d0,v1-v2); Lh/=np.linalg.norm(Lh); e1=d0/np.linalg.norm(d0); e2=np.cross(Lh,e1)
t=0.0; T=st['period']; tend=(orbits+0.5)*T; apo_t=[]; apo_ph=[]
hist=[]; vt=[]; vph=[]
while t<tend:
    d=x1-x2; r=math.sqrt(d@d); om=math.sqrt(K/(mu*r**3)); dt=min(cap,2*math.pi/(h*om))
    F1,F2=forces(x1,x2,v1,v2)
    p1=mom(v1,m1)+F1*(0.5*dt); p2=mom(v2,m2)+F2*(0.5*dt)
    w1,w2=vel(p1,m1),vel(p2,m2); y1=x1+w1*dt; y2=x2+w2*dt
    G1,G2=forces(y1,y2,w1,w2)
    c1,c2=vel(p1+G1*0.5*dt,m1),vel(p2+G2*0.5*dt,m2)
    G1,G2=forces(y1,y2,c1,c2)
    p1=p1+G1*0.5*dt; p2=p2+G2*0.5*dt
    nx1,nx2,nv1,nv2=y1,y2,vel(p1,m1),vel(p2,m2); nt=t+dt
    g0=(x1-x2)@(v1-v2); g1=(nx1-nx2)@(nv1-nv2)
    if g0>0 and g1<=0:
        fr=g0/(g0-g1); dd=(x1-x2)*(1-fr)+(nx1-nx2)*fr
        apo_t.append(t+fr*dt); apo_ph.append(math.atan2(dd@e2,dd@e1))
    if g0>0 and g1<=0 and hist:
        pts=[hist[-1],(t,x1-x2),(nt,nx1-nx2)]
        th=np.unwrap([math.atan2(p[1]@e2,p[1]@e1) for p in pts]); rho=[np.linalg.norm(p[1]) for p in pts]
        xx=th-th[1]; A=np.vstack([np.ones(3),xx,xx**2]).T; a0,b,cq=np.linalg.solve(A,rho); vt.append(t); vph.append(th[1]-b/(2*cq))
    hist=[(t,x1-x2)]
    x1,x2,v1,v2,t=nx1,nx2,nv1,nv2,nt
ph=np.unwrap(apo_ph); k=np.polyfit(apo_t,ph,1)[0]
kv=np.polyfit(vt,np.unwrap(vph),1)[0]; print('vertex Omega %.10e'%kv); print('replica h %g orbits %d cap %g: Omega %.10e  apo %d  first apo t %.12e'%(h,orbits,cap,k,len(apo_t),apo_t[0]))

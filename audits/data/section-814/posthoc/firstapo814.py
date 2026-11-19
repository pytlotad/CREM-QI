# POST HOC 814 (po wynikach): import sys, math, numpy as np
sys.path.insert(0,'.')
from pred814 import parse_start, darwin_on_first
from scipy.integrate import solve_ivp
f=sys.argv[1]; st=parse_start(f)
m1,m2,kq,c=st['m1'],st['m2'],st['kq'],st['c']
vel=lambda p,m: p/(m*math.sqrt(1+(p@p)/(m*m*c*c))); mom=lambda v,m: v*m/math.sqrt(1-(v@v)/(c*c))
def rhs(t,y):
    x1,x2,p1,p2=y[0:3],y[3:6],y[6:9],y[9:12]; v1,v2=vel(p1,m1),vel(p2,m2); d=x1-x2; r=np.linalg.norm(d)
    F1=d*(kq/r**3); F2=-F1
    return np.concatenate([v1,v2,F1+darwin_on_first(v1,v2,F2/m2,d,kq,c),F2+darwin_on_first(v2,v1,F1/m1,-d,kq,c)])
y0=np.concatenate([st['x1'],st['x2'],mom(st['v1'],m1),mom(st['v2'],m2)])
def ev(t,y): return (y[0:3]-y[3:6])@(vel(y[6:9],m1)-vel(y[9:12],m2))
ev.direction=-1
T=st['period']
sol=solve_ivp(rhs,(0,3.2*T),y0,method='DOP853',rtol=1e-12,atol=1e-30,events=ev,max_step=T/200)
for t,y in zip(sol.t_events[0],sol.y_events[0]): print('pred apo t %.12e  r %.9e %.9e  |r| %.12e'%(t,*(y[0:2]-y[3:5]),np.linalg.norm(y[0:3]-y[3:6])))
# engine
for ln in open(f):
    if ln.startswith('PAIR APSE808 mechB0 apo') and int(ln.split()[4])<=3:
        w=ln.split(); i=w.index('tc'); j=w.index('r1'); print('eng  apo t %s  r1 %s %s'%(w[i+1],w[j+1],w[j+2]))
# energies: Darwin-Lagrangian energy at start

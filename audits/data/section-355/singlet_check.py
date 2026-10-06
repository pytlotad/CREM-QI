# Audit 355: (1) QM -- the singlet is an eigenvector of S1.T.S2 for every SYMMETRIC T (control: non-symmetric T breaks it);
# (2) classical -- dS_i/dt = S_i x T S_j: antiparallel spins leave S1+S2 = 0 unless T is isotropic or s lies on a principal
# axis; parallel spins turn rigidly.  Pure Python (no numpy).
import random,math
random.seed(356)
def kron(a,b): return [[a[i//2][j//2]*b[i%2][j%2] for j in range(4)] for i in range(4)]
def mm(a,b): return [[sum(a[i][k]*b[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
def add(a,b,c=1): return [[a[i][j]+c*b[i][j] for j in range(4)] for i in range(4)]
def mv(a,v): return [sum(a[i][k]*v[k] for k in range(4)) for i in range(4)]
h=0.5; sx=[[0,h],[h,0]]; sy=[[0,-h*1j],[h*1j,0]]; sz=[[h,0],[0,-h]]; I=[[1,0],[0,1]]
S1=[kron(s,I) for s in (sx,sy,sz)]; S2=[kron(I,s) for s in (sx,sy,sz)]
sing=[0,1/math.sqrt(2),-1/math.sqrt(2),0]
def H_of(T):
    H=[[0]*4 for _ in range(4)]
    for i in range(3):
        for j in range(3): H=add(H,mm(S1[i],S2[j]),T[i][j])
    return H
def dev(H):
    v=mv(H,sing); e=sum(sing[i]*v[i] for i in range(4))
    return math.sqrt(sum(abs(v[i]-e*sing[i])**2 for i in range(4)))
ws=wa=0
for k in range(200):
    A=[[random.gauss(0,1) for _ in range(3)] for _ in range(3)]
    T=[[(A[i][j]+A[j][i])/2 for j in range(3)] for i in range(3)]
    ws=max(ws,dev(H_of(T))); wa=max(wa,dev(H_of(A)))
print(f"QM, symmetric T (200 random):   max |H s - <H> s| = {ws:.2e}")
print(f"QM, NON-symmetric T (control):  max |H s - <H> s| = {wa:.2e}")
def cross(a,b): return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def Tv(T,v): return [sum(T[i][j]*v[j] for j in range(3)) for i in range(3)]
def rhs(y,T): a,b=y[:3],y[3:]; return cross(a,Tv(T,b))+cross(b,Tv(T,a))
for label,T,s0 in [("isotropic T (contact)",[[1,0,0],[0,1,0],[0,0,1]],[1,0,1]),("tensor diag(1,1,-2), s at 45 deg",[[1,0,0],[0,1,0],[0,0,-2]],[1,0,1]),("tensor, s along principal axis",[[1,0,0],[0,1,0],[0,0,-2]],[0,0,1])]:
    for sign,name in [(-1,"antiparallel"),(1,"parallel")]:
        n=math.sqrt(sum(x*x for x in s0)); s=[x/n/2 for x in s0]; y=s+[sign*x for x in s]; dt=1e-3; mx=0; m0=math.sqrt(sum((y[i]+y[3+i])**2 for i in range(3)))
        for _ in range(20000):
            k1=rhs(y,T); k2=rhs([y[i]+dt/2*k1[i] for i in range(6)],T); k3=rhs([y[i]+dt/2*k2[i] for i in range(6)],T); k4=rhs([y[i]+dt*k3[i] for i in range(6)],T)
            y=[y[i]+dt/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]) for i in range(6)]
            mx=max(mx,abs(math.sqrt(sum((y[i]+y[3+i])**2 for i in range(3)))-m0))
        print(f"classical {label:34s} {name:12s}: max change of |S1+S2| = {mx:.3e}")

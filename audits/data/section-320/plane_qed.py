# o-Ps(xi) -> 3 gamma, tree-level QED, threshold. Rate as a function of the triplet
# polarization xi relative to the decay-plane normal n. Units m = 1.
import numpy as np, itertools
I2=np.eye(2); Z=np.zeros((2,2))
sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]]); sz=np.array([[1,0],[0,-1]],complex)
g0=np.block([[I2,Z],[Z,-I2]]).astype(complex)
gi=[np.block([[Z,s],[-s,Z]]) for s in (sx,sy,sz)]
def slash(a):  # a = (a0, ax, ay, az), metric (+---)
    return a[0]*g0 - a[1]*gi[0] - a[2]*gi[1] - a[3]*gi[2]
def mdot(a,b): return a[0]*b[0]-a[1:]@b[1:]
p=np.array([1.,0,0,0]); one=np.eye(4)
def S(q): return (slash(q)+one)/(mdot(q,q)-1.0)
def pols(k):
    kh=k/np.linalg.norm(k); t=np.array([0,0,1.]) if abs(kh[2])<0.9 else np.array([1.,0,0])
    e1=np.cross(kh,t); e1/=np.linalg.norm(e1); e2=np.cross(kh,e1); return [e1,e2]
def rate(ks,xi):
    Pi=slash(np.r_[0,xi])@(one-g0)
    tot=0.0
    K=[np.r_[np.linalg.norm(k),k] for k in ks]
    E=[pols(k) for k in ks]
    for pol in itertools.product(range(2),repeat=3):
        eps=[np.r_[0,E[i][pol[i]]] for i in range(3)]
        O=np.zeros((4,4),complex)
        for a,b,c in itertools.permutations(range(3)):
            O+=slash(eps[c])@S(p-K[a]-K[b])@slash(eps[b])@S(p-K[a])@slash(eps[a])
        tot+=abs(np.trace(O@Pi))**2
    return tot
rng=np.random.default_rng(1)
A=B=0.0; n=0; dev=[]
while n<4000:
    x1,x2=rng.uniform(0,1,2); x3=2-x1-x2
    if not (0<x3<1): continue
    # triangle with sides x1,x2,x3 in xy plane
    c12=(x3**2-x1**2-x2**2)/(2*x1*x2)
    k1=x1*np.array([1.,0,0]); k2=x2*np.array([c12,np.sqrt(max(0,1-c12**2)),0]); k3=-(k1+k2)
    ks=[k1,k2,k3]
    wz=rate(ks,np.array([0,0,1.])); wx=rate(ks,np.array([1.,0,0])); wy=rate(ks,np.array([0,1.,0]))
    op=((1-x1)/(x2*x3))**2+((1-x2)/(x1*x3))**2+((1-x3)/(x1*x2))**2
    dev.append((wx+wy+wz)/op)
    A+=0.5*(wx+wy); B+=wz; n+=1
dev=np.array(dev)
print("Ore-Powell check: spin-avg rate / OP formula: mean %.6g, rel spread %.2e"%(dev.mean(),dev.std()/dev.mean()))
print("pointwise wz/(wx+wy+wz) range: n/a; integrated a(in-plane)=%.6g b(normal)=%.6g b/a=%.6f"%(A/n,B/n,B/A))
r=B/A
# W(n|xi) = a(1-|xi.n|^2)+b|xi.n|^2.  m=+1 about s: |xi.n|^2=(1-c^2)/2 ; m=0: c^2
# m=+1: a(1+c^2)/2 + b(1-c^2)/2 ∝ (a+b) + (a-b)c^2
c2m1=(1-r)/(1+r); print("m=+-1: W ∝ 1 + %.6f cos^2 ; a2(P2) = %.6f"%(c2m1, (2/3)*c2m1/(1+c2m1/3)))
c20=(r-1); print("m=0  : W ∝ 1 + %.6f cos^2 ; a2(P2) = %.6f"%(c20/1, (2/3)*c20/(1+c20/3)))

# Audit 588: o-Ps -> 3 gamma plane law for ANY triplet polarization xi, tree-level QED (amplitude of audit 320a), complex tensor.
import numpy as np, itertools
src=open('../section-320/plane_qed.py').read().split('rng=np.random')[0]
g={}; exec(src,g)
slash,S,pols,one,g0,p=g['slash'],g['S'],g['pols'],g['one'],g['g0'],g['p']
def amps(ks):
    """M[pol, i]: amplitude for photon polarizations pol and triplet Cartesian component i (xi = e_i)."""
    K=[np.r_[np.linalg.norm(k),k] for k in ks]; E=[pols(k) for k in ks]
    out=[]
    for pol in itertools.product(range(2),repeat=3):
        eps=[np.r_[0,E[i][pol[i]]] for i in range(3)]
        O=np.zeros((4,4),complex)
        for a,b,c in itertools.permutations(range(3)):
            O+=slash(eps[c])@S(p-K[a]-K[b])@slash(eps[b])@S(p-K[a])@slash(eps[a])
        out.append([np.trace(O@slash(np.r_[0,e])@(one-g0)) for e in np.eye(3)])
    return np.array(out)
def rot(phi):
    c,s=np.cos(phi),np.sin(phi); return np.array([[c,-s,0],[s,c,0],[0,0,1.]])
rng=np.random.default_rng(588)
asym=[]; Tavg=np.zeros((3,3),complex); n=0
while n<400:
    x1,x2=rng.uniform(0,1,2); x3=2-x1-x2
    if not (0<x3<1): continue
    c12=(x3**2-x1**2-x2**2)/(2*x1*x2)
    k1=x1*np.array([1.,0,0]); k2=x2*np.array([c12,np.sqrt(max(0,1-c12**2)),0]); ks=[k1,k2,-(k1+k2)]
    M=amps(ks); T=M.T@M.conj()          # T_ij = sum_pol M_i M_j^*
    asym.append(np.max(np.abs(T-T.conj().T))/np.max(np.abs(T)) if False else np.max(np.abs(T.imag))/np.max(np.abs(T)))
    # in-plane average: rotate the triangle by phi_q (T -> R T R^T), 8 nodes
    for q in range(8):
        R=rot(2*np.pi*q/8); Tavg+=R@T@R.T/8
    n+=1
Tavg/=n
print(f"T1a  max_points |Im T|/|T| = {max(asym):.2e}")
A=0.5*(Tavg[0,0]+Tavg[1,1]).real; B=Tavg[2,2].real
off=np.max(np.abs(Tavg-np.diag(np.diag(Tavg))))/A
print(f"T1b  averaged T: A (in-plane) = {A:.6e}, B (normal) = {B:.6e}, B/A = {B/A:.12f}, |Txx-Tyy|/A = {abs(Tavg[0,0]-Tavg[1,1])/A:.1e}, max off-diag/A = {off:.1e}")
def W(xi,nh): return (xi.conj()@Tavg@xi).real
xis={'m=0 (real z)':np.array([0,0,1.],complex),'m=+1 (circular about z)':np.array([1,1j,0])/np.sqrt(2),'m=-1':np.array([1,-1j,0])/np.sqrt(2)}
for k in range(3):
    v=rng.normal(size=3)+1j*rng.normal(size=3); xis[f'elliptic {k}']=v/np.linalg.norm(v)
worst=0.0
print("T1c  W(n|xi)/[A(|xi|^2 + |xi.n|^2)] over 21 axis tilts (xi rotated, n = z fixed):")
for name,xi0 in xis.items():
    vals=[]
    for th in np.linspace(0,np.pi,21):
        Ry=np.array([[np.cos(th),0,np.sin(th)],[0,1,0],[-np.sin(th),0,np.cos(th)]])
        xi=Ry@xi0; pred=A*(np.vdot(xi,xi).real+abs(xi[2])**2)
        vals.append(W(xi,None)/pred)
    vals=np.array(vals); worst=max(worst,np.max(np.abs(vals-1)))
    print(f"     {name:24s} min {vals.min():.12f} max {vals.max():.12f}")
print(f"     max |ratio - 1| = {worst:.1e}")
c=np.linspace(-1,1,200001)
for lab,f in (("m=+-1: 1 + (1-c^2)/2",1+(1-c*c)/2),("m=0: 1 + c^2",1+c*c)):
    print(f"     <P2(n.axis)> for {lab}: {np.trapezoid(f*0.5*(3*c*c-1),c)/np.trapezoid(f,c):.6f}")

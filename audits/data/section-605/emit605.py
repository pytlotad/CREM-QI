# Audit 605: Born weights from emission energy balance + rule-Q single-moment projections (1S exact?), and why not for 2P.
import numpy as np, math
src=open('../section-601/qm601.py').read().split('print("B[T]')[0]; g={}; exec(src,g)
H2P,H1S,mu,h,dE,PT=g['H2P'],g['H1S'],g['mu'],g['h'],g['dE'],g['PT']
rng=np.random.default_rng(605)
def iso(n): v=rng.normal(size=(n,3)); return v/np.linalg.norm(v,axis=1,keepdims=True)
def spinor(nv):
    th=math.acos(max(-1,min(1,nv[2]))); ph=math.atan2(nv[1],nv[0]); return np.array([math.cos(th/2),np.exp(1j*ph)*math.sin(th/2)])
def lcoh(nv):  # l = 1 coherent state along nv, m = +1, 0, -1 (Wigner d^1_{m,1})
    c=nv[2]; ph=math.atan2(nv[1],nv[0]); s=math.sqrt(max(0,1-c*c))
    return np.array([np.exp(-1j*ph)*(1+c)/2, s/math.sqrt(2), np.exp(1j*ph)*(1-c)/2])
print("T1  1S: Born vs rule-Q projections + energy balance")
worst=0
for B in (0.01,0.03,0.1,0.3,1.0,3.0):
    H=H1S(B).real; w,V=np.linalg.eigh(H); K=dE/4; a=mu*B/h
    idxp=[i for i in range(4) if abs(V[0,i])**2>0.5][0]; idxm=[i for i in range(4) if abs(V[3,i])**2>0.5][0]; blk=[i for i in range(4) if i not in (idxp,idxm)]
    for n1,n2 in zip(iso(2000),iso(2000)):
        psi=np.kron(spinor(n1),spinor(n2)); born=np.abs(V.T.conj()@psi)**2
        pp=(1+n1[2])*(1+n2[2])/4; pm=(1-n1[2])*(1-n2[2])/4; P0=1-pp-pm
        Hcl=K*np.dot(n1,n2)+a*(n1[2]-n2[2])
        Erem=Hcl-pp*w[idxp]-pm*w[idxm]; Ea,Eb=w[blk[0]],w[blk[1]]
        pa=(Erem-P0*Eb)/(Ea-Eb); pb=P0-pa
        rec=np.zeros(4); rec[idxp]=pp; rec[idxm]=pm; rec[blk[0]]=pa; rec[blk[1]]=pb
        worst=max(worst,np.max(np.abs(rec-born)))
print(f"    max |reconstruction - Born| = {worst:.1e} (6 fields x 2000)")
print("T2  2P: block sizes and MaxEnt(P_M, <H>) vs Born")
Mop=None
def Mvals():
    # J_z eigenvalue of basis |m_L>|m1>|m2>, index 4 i + 2 p + q
    return np.array([ (1-i) + (0.5 if p==0 else -0.5) + (0.5 if q==0 else -0.5) for i in range(3) for p in range(2) for q in range(2)])
Mb=Mvals()
for B in (0.01,0.03,0.1,0.3,1.0,3.0):
    H=H2P(B); w,V=np.linalg.eigh(H); Mk=np.array([np.round(np.sum(np.abs(V[:,k])**2*Mb)) for k in range(12)])
    tri=np.array([(V[:,k].conj()@PT@V[:,k]).real for k in range(12)])
    sizes={int(m):int(np.sum(Mk==m)) for m in sorted(set(Mk))}
    tv=[]; d3=[]
    for Lh,n1,n2 in zip(iso(500),iso(500),iso(500)):
        psi=np.kron(np.kron(lcoh(Lh),spinor(n1)),spinor(n2)); born=np.abs(V.T.conj()@psi)**2
        PM={m:np.sum(born[Mk==m]) for m in sizes}                 # = rule-Q convolution of single-moment projections (checked below)
        Hmean=np.sum(born*w)
        def rec(beta):
            p=np.zeros(12)
            for m in sizes:
                sel=Mk==m; z=-beta*w[sel]; x=np.exp(z-z.max()); p[sel]=PM[m]*x/x.sum()
            return p
        lo,hi=-1e-8,1e-8
        for _ in range(200):
            if np.sum(rec(lo)*w)<Hmean: lo*=2
            else: break
        for _ in range(200):
            if np.sum(rec(hi)*w)>Hmean: hi*=2
            else: break
        for _ in range(100):
            mid=0.5*(lo+hi)
            if np.sum(rec(mid)*w)>Hmean: lo=mid
            else: hi=mid
        p=rec(0.5*(lo+hi)); tv.append(0.5*np.sum(np.abs(p-born))); d3.append(abs(np.sum(p*tri)-np.sum(born*tri)))
    print(f"    B {B:5.2f} T: block sizes {sizes}; total variation mean {np.mean(tv):.3f} max {np.max(tv):.3f}; |3P share MaxEnt - Born| mean {np.mean(d3):.3f} max {np.max(d3):.3f}")
# check P_M = convolution of single-moment rule-Q projections
Lh,n1,n2=iso(1)[0],iso(1)[0],iso(1)[0]
pl=np.abs(lcoh(Lh))**2; p1=np.array([(1+n1[2])/2,(1-n1[2])/2]); p2=np.array([(1+n2[2])/2,(1-n2[2])/2])
conv={}
for i,ml in enumerate((1,0,-1)):
    for a_,m1 in enumerate((0.5,-0.5)):
        for b_,m2 in enumerate((0.5,-0.5)): conv[ml+m1+m2]=conv.get(ml+m1+m2,0)+pl[i]*p1[a_]*p2[b_]
psi=np.kron(np.kron(lcoh(Lh),spinor(n1)),spinor(n2)); direct={}
for idx in range(12): direct[Mb[idx]]=direct.get(Mb[idx],0)+abs(psi[idx])**2
print("    P_M: rule-Q convolution vs |<M|psi>|^2 max diff", max(abs(conv[m]-direct[m]) for m in conv))
print("T3  classical energy vs <psi|H|psi>")
d1=0
for n1,n2 in zip(iso(500),iso(500)):
    B=0.3; K=dE/4; a=mu*B/h; psi=np.kron(spinor(n1),spinor(n2))
    d1=max(d1,abs((psi.conj()@H1S(B)@psi).real-(K*np.dot(n1,n2)+a*(n1[2]-n2[2])))/K)
print(f"    1S (spin 1/2, sigma = n): max relative diff {d1:.1e}")

# T3, 2P part: Q-symbol <L^, n1, n2|H|L^, n1, n2> vs classical H with |L| = 1 and 3/2, |s| = 1/2 (sigma = n), coefficients A, B from E_J
EJ=g['EJ']; A=(EJ[2]-EJ[1]+0.4*0)  # placeholder, replaced below
# from E_J = A <L.S> + B <S12>: solve with <L.S> = -2, -1, 1 and <S12> = -4, 2, -0.4 (least squares, exact)
M=np.array([[-2,-4],[-1,2],[1,-0.4]]); A,Bt=np.linalg.lstsq(M,EJ,rcond=None)[0]
for B in (0.1,1.0):
    H=H2P(B); a=mu*B/h; dev={1.0:0,1.5:0}; scale=np.max(np.abs(np.linalg.eigvalsh(H)))
    for Lh,n1,n2 in zip(iso(300),iso(300),iso(300)):
        psi=np.kron(np.kron(lcoh(Lh),spinor(n1)),spinor(n2)); q=(psi.conj()@H@psi).real
        for Llen in (1.0,1.5):
            L=Llen*Lh; s1=0.5*n1; s2=0.5*n2; p1=np.dot(s1,Lh); p2=np.dot(s2,Lh); d12=np.dot(s1,s2)
            Hcl=A*np.dot(L,s1+s2)+4*Bt*(1.5*(d12-p1*p2)-d12)+2*a*(s1[2]-s2[2])
            dev[Llen]=max(dev[Llen],abs(q-Hcl)/scale)
    print(f"    2P B={B} T: max |<psi|H|psi> - H_cl|/max|E|: |L| = 1: {dev[1.0]:.3f}, |L| = 3/2: {dev[1.5]:.3f}  (A = {A:.4f}, B = {Bt:.5f} GHz)")

# Audit 492 A: author's family sigma_1,2 = a S^ +- u (491f-j) versus the engine (sigma = 1, triplet u = 0).
# Planar orbit r = (cos p, sin p, 0), L || z; S^ = (sin b, 0, cos b); u latched at alpha from the projection of r onto the
# plane perpendicular to S^ (as in s12_zatrzask.py); "fixed u" = u fixed in space (no latch), orbit-averaged.
# Swap the ingredients one at a time from the author's construction to the engine.  python3 s492.py (numpy).
import numpy as np
P=np.linspace(0,2*np.pi,1441)[:-1]
R=np.stack([np.cos(P),np.sin(P),0*P],1)
def s12(sig,a,c,alpha,latched=True):
    """<3(s1.r)(s2.r) - s1.s2> over the orbit; |s_i| = sig, projection a on S^ at cos c to L."""
    u2=sig*sig-a*a
    if u2<-1e-12: return np.nan
    u=np.sqrt(max(u2,0.0)); s=np.sqrt(max(0.0,1-c*c)); Sh=np.array([s,0,c])
    e1=np.array([c,0,-s]); e2=np.cross(Sh,e1)
    if latched:
        rp=R-np.outer(R@Sh,Sh); n=np.linalg.norm(rp,axis=1)
        d=np.where(n[:,None]>1e-12,rp/np.maximum(n,1e-300)[:,None],e1); d2=np.cross(Sh,d)
        U=u*(np.cos(alpha)*d+np.sin(alpha)*d2)
    else:
        U=np.tile(u*(np.cos(alpha)*e1+np.sin(alpha)*e2),(len(P),1))
    s1=a*Sh+U; s2=a*Sh-U
    return np.mean(3*np.sum(s1*R,1)*np.sum(s2*R,1)-np.sum(s1*s2,1))
AL=np.radians(np.arange(0,360,2.0))
def rng(sig,a,c,latched=True):
    v=[s12(sig,a,c,x,latched) for x in AL]; return min(v),max(v),np.mean(v)
def tri(L,S,J): return (J*J-L*L-S*S)/(2*L*S)
QM={0:-4,1:2,2:-0.4}
lenQM=dict(L=np.sqrt(2),S=np.sqrt(2),sig=np.sqrt(3),J=lambda j:np.sqrt(j*(j+1)))
lenLa=dict(L=1.5,S=1.5,sig=2.0,J=lambda j:j+0.5)
print("A1. ingredient swaps per J (cos L-S from the triangle of the named lengths); [min, max] over alpha, mean over alpha")
rows=[("author QM: sig sqrt3, a sqrt2, latch", lenQM, lambda d:(d['sig'],d['S'],True)),
      ("author Langer: sig 2, a 1.5, latch",   lenLa, lambda d:(d['sig'],d['S'],True)),
      ("(iv) QM lengths, u fixed in space",    lenQM, lambda d:(d['sig'],d['S'],False)),
      ("(ii) QM: sig = a = sqrt2 (u = 0)",     lenQM, lambda d:(d['S'],d['S'],True)),
      ("(i)  sig 1, a 1 (u = 0), QM angle",    lenQM, lambda d:(1.0,1.0,True)),
      ("(i)  sig 1, a 1 (u = 0), Q angle",     lenLa, lambda d:(1.0,1.0,True)),
      ("     sig 1, a 0.5, latch, Q angle",    lenLa, lambda d:(1.0,0.5,True))]
for name,d,f in rows:
    print("  "+name)
    for J in (0,1,2):
        c=np.clip(tri(d['L'],d['S'],d['J'](J)),-1,1); sig,a,lat=f(d); lo,hi,m=rng(sig,a,c,lat)
        hit='yes' if lo-1e-9<=QM[J]<=hi+1e-9 else 'no'
        print(f"    J={J} cos(L,S) {c:+.3f}  [{lo:+.3f}, {hi:+.3f}]  mean {m:+.3f}   QM {QM[J]:+}  reachable {hit}")
print("A2. global range over cos(L,S) in [-1,1] (step 0.02), alpha (2 deg) and, for sig = 1, a in [0,1] (step 0.05)")
C=np.linspace(-1,1,101)
def glob(sig,As,lat=True):
    best=(9,None); worst=(-9,None)
    for a in As:
        for c in C:
            for x in AL[::2]:
                v=s12(sig,a,c,x,lat)
                if v<best[0]: best=(v,(a,c,np.degrees(x)))
                if v>worst[0]: worst=(v,(a,c,np.degrees(x)))
    return best,worst
for name,sig,As,lat in [("author QM sig sqrt3, a sqrt2, latched",np.sqrt(3),[np.sqrt(2)],True),
                        ("author Langer sig 2, a 1.5, latched",2.0,[1.5],True),
                        ("engine triplet sig 1, a 1 (u = 0)",1.0,[1.0],True),
                        ("engine free spins sig 1, a in [0,1], latched",1.0,np.linspace(0,1,21),True),
                        ("engine free spins sig 1, a in [0,1], u fixed",1.0,np.linspace(0,1,21),False)]:
    (lo,al),(hi,ah)=glob(sig,As,lat)
    print(f"  {name}: min {lo:+.3f} at (a, cos, alpha) = ({al[0]:.2f}, {al[1]:+.2f}, {al[2]:.0f})  max {hi:+.3f} at ({ah[0]:.2f}, {ah[1]:+.2f}, {ah[2]:.0f})")
print("A3. engine's own classical tensor (sig 1, a 1, u = 0) averaged over its rule-Q rounding basins, isotropic S^ (uniform cos)")
cs=np.linspace(-1,1,20001); Jl=np.sqrt(1.5**2+1.5**2+2*1.5*1.5*cs); Jr=np.clip(np.rint(Jl-0.5),0,2)
val=1.5*(1-cs**2)-1
for J in (0,1,2):
    m=Jr==J; print(f"  J={J}: basin cos in [{cs[m].min():+.3f}, {cs[m].max():+.3f}], weight {m.mean():.3f} (QM {(2*J+1)/9:.3f}), <S12> {val[m].mean():+.3f}  (QM {QM[J]:+})")
print("A4. identity check at J = 0 (S^ = -z, u || r): a^2 - 2 sig^2")
for sig,a in [(np.sqrt(3),np.sqrt(2)),(2.0,1.5),(1.0,1.0),(1.0,0.0)]:
    print(f"  sig {sig:.3f} a {a:.3f}: numeric {s12(sig,a,-1.0,0.0):+.4f}  formula {a*a-2*sig*sig:+.4f}")

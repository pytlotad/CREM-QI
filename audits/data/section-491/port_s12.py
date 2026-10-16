# literal pure-Python port of s12_zatrzask.py (same grids: 720 orbit azimuths, 91 latch angles on [0, pi/2], 72 for the unlatched mean)
import math
def dot(a,b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def cross(a,b): return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def add(a,b): return [a[i]+b[i] for i in range(3)]
def sub(a,b): return [a[i]-b[i] for i in range(3)]
def mul(a,k): return [x*k for x in a]
def norm(a): return math.sqrt(dot(a,a))
def run(Ll,Sl,sig,Jl):
    x=(Jl**2-Ll**2-Sl**2)/2; c=x/(Ll*Sl); s=math.sqrt(max(0,1-c*c))
    Sh=[s,0,c]; a=Sl
    u2=sig**2-a**2
    if u2<0: return None
    e1=[c,0,-s]; e2=cross(Sh,e1)
    phis=[2*math.pi*i/720 for i in range(720)]
    def S12(alpha):
        tot=0
        for p in phis:
            r=[math.cos(p),math.sin(p),0]
            rp=sub(r,mul(Sh,dot(r,Sh))); n=norm(rp)
            d=e1 if n<1e-12 else mul(rp,1/n)
            d2=cross(Sh,d)
            u=mul(add(mul(d,math.cos(alpha)),mul(d2,math.sin(alpha))),math.sqrt(u2))
            s1=add(mul(Sh,a),u); s2=sub(mul(Sh,a),u)
            tot+=3*dot(s1,r)*dot(s2,r)-dot(s1,s2)
        return tot/len(phis)
    al=[ (math.pi/2)*i/90 for i in range(91)]; v=[S12(x_) for x_ in al]
    vav=sum(S12(2*math.pi*i/72) for i in range(72))/72
    return x,v[0],v[-1],vav,(al,v)
QM={0:-4,1:2,2:-0.4}
for name,Lf,Sf,sigf,Jf in [("QM lengths",lambda l:math.sqrt(2),lambda:math.sqrt(2),math.sqrt(3),lambda J:math.sqrt(J*(J+1))),
                          ("Langer/rule Q",lambda l:1.5,lambda:1.5,2.0,lambda J:J+0.5)]:
    print(name)
    for J in (0,1,2):
        x,v0,v90,vav,(al,v)=run(Lf(1),Sf(),sigf,Jf(J))
        q=QM[J]
        ang=[math.degrees(al[i]) for i in range(len(v)-1) if (v[i]-q)*(v[i+1]-q)<0 or v[i]==q]
        print(f"  J={J} L.S={x:.3f}  latch u||r: {v0:.3f}  u_|_r: {v90:.3f}  unlatched {vav:.3f}  QM {q}  latch angle hitting QM: {[round(a_,1) for a_ in ang]}  range [{min(v):.3f}, {max(v):.3f}]")

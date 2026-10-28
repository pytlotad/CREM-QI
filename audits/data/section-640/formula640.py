"""Audit 640: Ps n = 2 levels from Pineda & Soto, PRD 59, 016005 (1998), hep-ph/9805424, eq. (3.8), against the engine (out/p640_*.txt) and [K04].
ln R = ln k0 (Bethe logarithms, hydrogen values); checked first on 1^3S1 - 2^3S1 alpha^5 = -1501.44 MHz = Karshenboim IJMPA 19, 3879 (2004) table 5."""
import math
a=7.2973525693e-3; mc2=1.2355899662e20/1e6; P=mc2*a**4/8; ap=a/math.pi
lnk0={(1,0):2.984128556,(2,0):2.811769893,(2,1):-0.030016709}
H=lambda n: sum(1/k for k in range(1,n+1))
def SI(n,l):
    if l==0: return P*(-(2*a/(3*math.pi))/n**3*(9*math.log(a)+7*math.log(n)+8*lnk0[(n,0)]-14*math.log(2)-49/15-7*(H(n)+(n-1)/(2*n))))
    return P*(-(16*a/(3*math.pi))/n**3*(lnk0[(n,l)]+7/16/(l*(l+1)*(2*l+1))))
spinS=lambda n: P*(14/3)/n**3*(3*a/(7*math.pi))*(-32/9-2*math.log(2))
print("check 1^3S1-2^3S1 alpha^5: %.2f MHz (K04 -1501.4)"%(SI(2,0)+spinS(2)-SI(1,0)-spinS(1)))
n=2; base=lambda l: P*(-4/(n**3*(2*l+1))+11/(8*n**4))
S4=base(0)+P*(14/3)/n**3; C0={0:-8,1:-2,2:2.8}; C1={0:-6*ap,1:-ap,2:1.8*ap}
P4={J:base(1)+P*C0[J]/(6*n**3) for J in range(3)}; P5={J:P*C1[J]/(6*n**3) for J in range(3)}
def read(f):
    d={}
    for line in open(f):
        w=line.split(); d[w[0]]=float(w[1])
    return d
e0=read("out/p640_V0.txt"); e1=read("out/p640_V1_g2_LO.txt")
print("alpha^4 formula - engine(g=2, LO): S %.4f P %s 1P1 %.4f"%(S4-e1["2^3S1"],[round(P4[J]-e1["2^3P%d"%J],4) for J in range(3)],base(1)-e1["2^1P1"]))
print("alpha^5 spin: S formula %.4f engine %.4f; P formula %s engine %s"%(spinS(2),e0["2^3S1"]-e1["2^3S1"],[round(P5[J],4) for J in range(3)],[round(e0["2^3P%d"%J]-e1["2^3P%d"%J],4) for J in range(3)]))
print("alpha^5 spin-independent (not in engine): 2S %.2f 2P %.2f MHz"%(SI(2,0),SI(2,1)))
K={"2^3P0":18498.25,"2^3P1":13012.41,"2^3P2":8625.70,"2^1P1":11185.37}
fS=S4+spinS(2)+SI(2,0)
for k,v in K.items():
    fP=(P4[int(k[-1])]+P5[int(k[-1])]+SI(2,1)) if k!="2^1P1" else base(1)+SI(2,1)
    eng=e0["2^3S1"]-e0[k]
    print("%s: engine %.2f formula %.2f K04 %.2f | engine-K04 %.2f engine-formula %.2f formula-K04 %.2f"%(k,eng,fS-fP,v,eng-v,eng-(fS-fP),(fS-fP)-v))

# Wagi skladowych Zeemana (import 20, audyt 505): klasyczne wagi Lorentza wobec Clebscha-Gordana
# dla przejscia E1 (l, m) -> (l-1, m'), z m' = m - q, q = +1 (sigma+), 0 (pi), -1 (sigma-).
# Lorentz: sigma (1 +- c)^2/4, pi (1 - c^2)/2, c = M/L.
# Warianty: (a) c z poczatkowego stanu przy L = l + 1/2 (Langer)
#           (b) srednia geometryczna wag klasycznych przy M poczatkowym i koncowym, L = l (srodek miedzy l+1/2 i l-1/2)
from fractions import Fraction as F
def qm(l,m,q):
    d=l*(2*l+1)*2
    if q==0: v=F(2*(l*l-m*m),d)
    elif q==1: v=F((l+m)*(l+m-1),d)
    else: v=F((l-m)*(l-m-1),d)
    return v/F(2*l-1,2*l+1)
def lor(c,q): return (1+c)**2/4 if q==1 else ((1-c)**2/4 if q==-1 else (1-c*c)/2)
def var_a(l,m,q): L=l+0.5; return lor(m/L,q)
def var_b(l,m,q):
    L=l; mf=m-q
    if q==0: return (1-(m/L)**2)/2   # M bez zmiany
    s=1 if q==1 else -1
    return ((1+s*m/L)*(1+s*mf/L))/4
worst={'a':0,'b':0}
for l in range(1,7):
    for m in range(-l,l+1):
        qs=[q for q in (1,0,-1) if abs(m-q)<=l-1]
        for key,f in (('a',var_a),('b',var_b)):
            raw={q:f(l,m,q) for q in qs}; s=sum(raw.values())
            for q in qs: worst[key]=max(worst[key],abs(raw[q]/s-float(qm(l,m,q))))
        if l==2 and m==0:
            print("3d m=0 -> 2p: QM",[round(float(qm(l,m,q)),4) for q in qs],
                  " Lorentz(a)",[round(var_a(l,m,q)/sum(var_a(l,m,x) for x in qs),4) for q in qs],
                  " srednia geom.(b)",[round(var_b(l,m,q)/sum(var_b(l,m,x) for x in qs),4) for q in qs])
print("najwieksze odchylenie od CG^2 dla l=1..6, wszystkie m:",{k:round(v,6) for k,v in worst.items()})

# Zasieg <S12> w 3P_J przy zatrzasku przeciwfazowej skladowej spinow na azymucie polozenia.
# Plaska orbita, stale dlugosci |L|,|S|,|sigma_i|, zachowane J. Dwa zestawy dlugosci: QM i Langer/regula Q.
# Wynik: tabela w sprzezenie-fazowe-l.md, sekcja 2. Uruchom: python3 s12_zatrzask.py (wymaga numpy).
import numpy as np
def run(Ll,Sl,sig,Jl,label):
    # lengths in hbar: |L|=Ll, |S|=Sl, |J|=Jl ; |sigma_i|=sig (sigma units, sigma=2 s/hbar)
    x=(Jl**2-Ll**2-Sl**2)/2; c=x/(Ll*Sl); s=np.sqrt(max(0,1-c*c))
    Sh=np.array([s,0,c]); a=Sl  # sigma1+sigma2 = 2*S/hbar -> each projection a=Sl
    u2=sig**2-a**2
    if u2<0: return None
    e1=np.array([c,0,-s]); e2=np.cross(Sh,e1)
    phis=np.linspace(0,2*np.pi,721)[:-1]
    def S12(alpha):  # u latched: angle alpha from projection of r onto plane perp S
        tot=0
        for p in phis:
            r=np.array([np.cos(p),np.sin(p),0])
            rp=r-(r@Sh)*Sh; n=np.linalg.norm(rp)
            if n<1e-12: d=e1
            else: d=rp/n
            d2=np.cross(Sh,d)
            u=np.sqrt(u2)*(np.cos(alpha)*d+np.sin(alpha)*d2)
            s1=a*Sh+u; s2=a*Sh-u
            tot+=3*(s1@r)*(s2@r)-s1@s2
        return tot/len(phis)
    al=np.linspace(0,np.pi/2,91); v=[S12(x_) for x_ in al]
    return x, v[0], v[-1], np.mean([S12(t) for t in np.linspace(0,2*np.pi,73)[:-1]]), (al,v)
QM={0:-4,1:2,2:-0.4}
for name,Lf,Sf,sigf,Jf in [("QM lengths",lambda l:np.sqrt(2),lambda:np.sqrt(2),np.sqrt(3),lambda J:np.sqrt(J*(J+1))),
                          ("Langer/rule Q",lambda l:1.5,lambda:1.5,2.0,lambda J:J+0.5)]:
    print(name)
    for J in (0,1,2):
        x,v0,v90,vav,(al,v)=run(Lf(1),Sf(),sigf,Jf(J),name)
        q=QM[J]; vv=np.array(v)
        idx=np.where(np.diff(np.sign(vv-q)))[0]
        ang=[np.degrees(al[i]) for i in idx]
        print(f"  J={J} L.S={x:.3f}  latch u||r: {v0:.3f}  u_|_r: {v90:.3f}  unlatched {vav:.3f}  QM {q}  latch angle hitting QM: {ang}")

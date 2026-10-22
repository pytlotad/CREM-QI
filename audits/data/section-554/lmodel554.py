# Audit 554: model l' shares within n from FORCEDEVENT lines (E, L of the captured state), weights HT x P / sigma_Kramers(K) (as ana551).  python3 lmodel554.py stdout stderr B_eV pair
import sys, re, math, numpy as np
out,err,B,pair=sys.argv[1],sys.argv[2],float(sys.argv[3]),sys.argv[4]
eps0=8.8541878128e-12; e=1.602176634e-19; hbar=1.054571817e-34; me=9.1093837015e-31; mp_=1.67262192369e-27; alpha=7.2973525693e-3
k=e*e/(4*np.pi*eps0); mu=me/2 if pair=='Ps' else me*mp_/(me+mp_); a=hbar**2/(mu*k)
def sigb(K): Kj=K*e; rK=(hbar**2*k/(mu*Kj*Kj))**(1/3); return 0.5*math.sqrt(rK*(rK+k/Kj))
def kram(K):
    eta2=B/K; ns=np.arange(1,20000); return ((32*np.pi/(3*np.sqrt(3)))*alpha**3*a*a*eta2**2/(ns*(eta2+ns*ns))).sum()
W={}; tot=0; Wn={}
for l in open(err):
    if not l.startswith('FORCEDEVENT'): continue
    d=dict(x.split('=',1) for x in l.split()[1:])
    if d['captured'] not in ('Para-Positronium','Ortho-Positronium'): continue
    K=float(d['K']); b=float(d['b']); P=float(d['P']); E=float(d['E']); L=float(d['L'])
    if not (E<0): continue
    n=max(1,round(math.sqrt(B/abs(E)))); lp=min(n-1,max(0,round(L-0.5)))
    w=P*2*math.pi*sigb(K)**2*math.exp(b*b/(2*sigb(K)**2))/kram(K)
    W[(n,lp)]=W.get((n,lp),0)+w; Wn[n]=Wn.get(n,0)+w; tot+=w
for n in range(1,5):
    if n not in Wn: continue
    print(f"  n = {n}: share {Wn[n]/tot:.3f}; l' shares "+' '.join(f"l={l}: {W.get((n,l),0)/Wn[n]:.3f}" for l in range(n)))

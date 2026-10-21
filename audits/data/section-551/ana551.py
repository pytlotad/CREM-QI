# Audit 551 (from 550, energy-scaled sigma_b(K) = b_max(K)/2 per event): cross section and level shares from the forced-capture lines (FORCEDEVENT, stderr) and the CAPDUMP lines (stdout); Lambda against expect549's semi-analytic
# value event by event.  python3 ana550.py stdout stderr sigma_b_pm B_eV mu(Ps|H)
import sys, re, math
sys.path.insert(0,'../section-549')
import numpy as np
out,err,sb,B,pair=sys.argv[1],sys.argv[2],float(sys.argv[3])*1e-12,float(sys.argv[4]),sys.argv[5]
eps0=8.8541878128e-12; e=1.602176634e-19; hbar=1.054571817e-34; c=299792458.0; me=9.1093837015e-31; mp=1.67262192369e-27
k=e*e/(4*np.pi*eps0); mu=me/2 if pair=='Ps' else me*mp/(me+mp); a=hbar**2/(mu*k); Bj=k/(2*a); alpha=7.2973525693e-3
def N_allowed(b,K):
    v=np.sqrt(2*K/mu); ah=k/(2*K); T=np.sqrt(mu*ah**3/k); rlo=(hbar**2*k/(mu*(K+4*Bj)**2))**(1/3); rhi=(hbar**2*k/(mu*K*K))**(1/3)
    L=mu*v*b; ecc=np.sqrt(1+2*K*L*L/(mu*k*k)); rp=ah*(ecc-1)
    if rp>=rhi: return 0.0
    Fm=np.arccosh((rhi/ah+1)/ecc); F=np.linspace(-Fm,Fm,4001); u=ecc*np.cosh(F)-1; r=ah*u
    om=np.sqrt(k/(mu*r**3)); P=e*e*(k/(mu*r*r))**2/(6*np.pi*eps0*c**3); dN=P/(hbar*om)*T*u
    return np.trapezoid(dN*(r>=rlo),F)
ev=[]; mainbound=0
for l in open(out):
    if l.startswith('CAPDUMP '):
        o=re.search(r'out=(\S+)',l).group(1); K=float(re.search(r' K=(\S+)',l).group(1)); b=float(re.search(r' b=(\S+)',l).group(1))
        ev.append((o,K,b)); mainbound+=o in ('Para-Positronium','Ortho-Positronium')
forced={}
for l in open(err):
    if l.startswith('FORCEDEVENT'):
        d=dict(x.split('=',1) for x in l.split()[1:]); forced[(float(d['K']),float(d['b']))]=d
def sigb(K):
    Kj=K*e; rK=(hbar**2*k/(mu*Kj*Kj))**(1/3); return 0.5*math.sqrt(rK*(rK+k/Kj))
HT=lambda b,K: 2*math.pi*sigb(K)**2*math.exp(b*b/(2*sigb(K)**2))
classified=[x for x in ev if x[0] not in ('Unresolved','NumericalFailure')]
vals=[]; nW={}; ratio=[]; Ks=[]
for o,K,b in classified:
    d=forced.get((K,b)); v=0.0
    if d:
        P=float(d['P']); cap=d['captured'] in ('Para-Positronium','Ortho-Positronium'); E=float(d['E'])
        if cap and E<0:
            v=P; n=max(1,round(math.sqrt(B/abs(E)))); nW[n]=nW.get(n,0)+P*HT(b,K)
        Na=N_allowed(b,K*e)
        if Na>0: ratio.append(float(d['Lambda'])/Na)
    vals.append(v*HT(b,K)); Ks.append(K)
vals=np.array(vals); sig=vals.mean(); se=vals.std()/math.sqrt(len(vals))
Km=np.mean(Ks); eta2=B/Km; ns=np.arange(1,20000); sk=(32*np.pi/(3*np.sqrt(3)))*alpha**3*a*a*eta2**2/(ns*(eta2+ns*ns))
tot=sum(nW.values())
def kram(KeV):
    eta2=B/KeV; ns=np.arange(1,20000); return ((32*np.pi/(3*np.sqrt(3)))*alpha**3*a*a*eta2**2/(ns*(eta2+ns*ns))).sum()
kmean=np.mean([kram(K) for K in Ks])
print(f"{pair}: classified {len(classified)}, forced events {len(forced)}, photon-free-pass captures (residual exchange, excluded) {mainbound}")
print(f"  sigma (forced only) {sig/1e-28:.4e} +- {se/1e-28:.2e} b; Kramers at <K> = {Km:.3f} eV: {sk.sum()/1e-28:.4e} b; ratio {sig/sk.sum():.3f} +- {se/sk.sum():.3f}")
print(f"  Kramers averaged over the sampled energies {kmean/1e-28:.4e} b; ratio {sig/kmean:.3f} +- {se/kmean:.3f}")
print("  level shares (HT x P): "+' '.join(f"n={n}: {nW[n]/tot:.3f}" for n in sorted(nW) if n<=6)+f"  (n>6: {sum(v for n,v in nW.items() if n>6)/tot:.3f}) | Kramers: "+' '.join(f"{x:.3f}" for x in sk[:6]/sk.sum()))
r=np.array(ratio); print(f"  Lambda(engine)/N_allowed(hyperbola): median {np.median(r):.4f}, 10/90 % {np.percentile(r,10):.4f}/{np.percentile(r,90):.4f} over {len(r)} events")
# Audit 551: the ratio estimated event by event, <HT x P / sigma_Kramers(K_i)>, which removes the 1/K spread of the energy sample, and level shares with the same weights
# against the event-averaged Kramers shares.
rat=[]; nR={}; kshare=np.zeros(6); cnt=0
for o,K,b in classified:
    sK=kram(K); d=forced.get((K,b)); v=0.0
    eta2=B/K; ns=np.arange(1,20000); sk=(32*np.pi/(3*np.sqrt(3)))*alpha**3*a*a*eta2**2/(ns*(eta2+ns*ns)); kshare+=sk[:6]/sk.sum(); cnt+=1
    if d:
        P=float(d['P']); E=float(d['E'])
        if d['captured'] in ('Para-Positronium','Ortho-Positronium') and E<0:
            v=P*HT(b,K)/sK; n=max(1,round(math.sqrt(B/abs(E)))); nR[n]=nR.get(n,0)+v
    rat.append(v)
rat=np.array(rat); tR=sum(nR.values())
print(f"  <HT x P / sigma_Kramers(K)> = {rat.mean():.3f} +- {rat.std()/math.sqrt(len(rat)):.3f} (rule on the hyperbola: 0.84-0.95)")
print("  level shares (HT x P / sigma_K): "+' '.join(f"n={n}: {nR[n]/tR:.3f}" for n in sorted(nR) if n<=6)+f"  (n>6: {sum(v for n,v in nR.items() if n>6)/tR:.3f}) | Kramers (event mean): "+' '.join(f"{x:.3f}" for x in kshare/cnt))

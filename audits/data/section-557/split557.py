# Audit 557: the exact QM cascade started from the MODEL's capture distribution (samples, l = round(L/hbar - 1/2) clipped to [0, n-1]) -- separates capture errors from cascade errors.
import sys, os, math
src=open('cascade557.py').read().split("# capture shares at K/B")[0].replace('print("lifetimes [ns]: "','#')
ns={}; exec(src,ns); A=ns['A']; tot=ns['tot']
def cascade(init):
    flux=dict(init); lines={}
    for n in range(8,1,-1):
        for l in range(n):
            f=flux.get((n,l),0)
            if f==0 or (n,l)==(2,0): continue
            for (a,b,m,lp),v in A.items():
                if (a,b)!=(n,l): continue
                p=f*v/tot[(n,l)]; flux[(m,lp)]=flux.get((m,lp),0)+p; lines[(n,m)]=lines.get((n,m),0)+p
    return flux,lines
for samples,d,label in (('samples-r553.txt','r553','553'),('samples-r554.txt','r554','554')):
    S=[l.split() for l in open(samples)]; init={}; N=0
    base='/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r542/casc/'+d
    for i,n,L in S:
        n=int(n); L=float(L)
        if n>1:
            txt=open(os.path.join(base,f'c{i}.err')).read()
            if 'stop 0 ' in txt: continue        # same censored samples excluded as in ana557
        l=min(n-1,max(0,round(L-0.5))); init[(n,l)]=init.get((n,l),0)+1; N+=1
    init={k:v/N for k,v in init.items()}
    flux,lines=cascade(init)
    lya=lines.get((2,1),0); direct=sum(v for (a,m),v in lines.items() if m==1 and a>2)
    print(f"QM cascade from model capture {label} ({N} samples): direct 1s {init.get((1,0),0):.3f}; Lyman-alpha {lya:.3f}; ends in 2s {flux.get((2,0),0):.3f}; Lyman-beta+ {direct:.3f}; "
          f"Halpha {lines.get((3,2),0):.3f}, Hbeta {lines.get((4,2),0):.3f}, Halpha/Hbeta {lines.get((3,2),0)/max(lines.get((4,2),0),1e-300):.2f}")

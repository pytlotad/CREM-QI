# Audit 613: engine today vs A (today's total, G shares) vs B (all from G) -- branching vs QM, lifetimes vs measurements (hydrogen).
import numpy as np
src=open('../section-560/tempo_e1_punkt_srodkowy.py').read().split('if __name__ == "__main__":')[0]; g={}; exec(src,g)
rate_qm,rate_rule,C,AU,channels=g['rate_qm'],g['rate_rule'],g['C'],g['AU_RATE'],g['channels']
s599=open('../section-599/var599.py').read().split('rows=[]')[0]; g2={}; exec(s599,g2); rG=lambda *a: g2['rate']('G','o',*a)
def harmonics(n,l,N=1<<20):
    a=n*n; L=l+0.5; e=np.sqrt(max(0.0,1-(L/n)**2)); b=a*np.sqrt(1-e*e)
    M=np.arange(N)*2*np.pi/N; u=np.full_like(M,np.pi)
    for _ in range(200): u-=(u-e*np.sin(u)-M)/(1-e*np.cos(u))
    z=a*(np.cos(u)-e)+1j*b*np.sin(u); c=np.fft.fft(z)/N; k=np.arange(1,N//2); Om=1.0/n**3
    return k, 2/3*np.abs(c[k])**2*(k*Om)**3/C**3
def rates(n,l):
    ch=list(channels(n,l)); qm={c:rate_qm(n,l,*c) for c in ch}
    down=[c for c in ch if c[1]==l-1]; up=[c for c in ch if c[1]==l+1]
    k,r=harmonics(n,l); today={}
    for c in down:
        n2=c[0]; today[c]=r[(np.maximum(n-k,l))==n2].sum()
    for c in up: today[c]=rate_rule(n,l,*c,'RBw')
    G={c:rG(n,l,*c) for c in ch}
    totd=sum(today[c] for c in down); gd=sum(G[c] for c in down)
    A={c:(totd*G[c]/gd if c in down else today[c]) for c in ch}
    return qm,{'today':today,'A':A,'B':G}
meas={(2,1):(1.600,0.004),(3,1):(5.58,0.13),(4,1):(11.25,0.78),(5,1):(21.9,3.0)}
worst={v:(0,None) for v in ('today','A','B')}; chi={v:0.0 for v in worst}; rows=[]
for l in (1,2,3):
    for n in range(l+1,9):
        qm,var=rates(n,l); tq=sum(qm.values()); tauq=1e9/(tq*AU)
        line=f"{n}{'spdf'[l]}: tau_QM {tauq:8.3f}"
        for v,R in var.items():
            tv=sum(R.values()); tau=1e9/(tv*AU)
            d=max(abs(R[c]/tv-qm[c]/tq) for c in qm)
            if d>worst[v][0]: worst[v]=(d,f"{n}{'spdf'[l]}")
            line+=f"  {v}: tau/QM {tau/tauq:.3f} maxdev {d:.3f}"
            if (n,l) in meas: m,s=meas[(n,l)]; chi[v]+=((tau-m)/s)**2; line+=f" (tau {tau:.3f} vs {m}({s}))"
        print(line,flush=True)
print("max branching deviation:",{v:f"{worst[v][0]:.3f} at {worst[v][1]}" for v in worst})
print("chi^2 vs measurements (2p, 3p Bickel-Goodman; 4p, 5p Etherton):",{v:round(chi[v],2) for v in chi})

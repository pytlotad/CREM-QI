# Per-photon path for each trajectory: (n, L/hbar, e^2) before -> after, trimmed?
import glob,re,math,statistics as st
S=''  # run in audits/data/section-345
fam={}
for l in open(S+'families.txt'):
    k,*ids=l.split(); fam.update({int(i):k for i in ids})
recs=[]; cur=None
import os
for f in ['trace.txt']:
    for l in open(S+f):
        if l.startswith('BEGIN'): _,i,ph=l.split(); cur={'i':int(i),'ch':'p-Ps' if ph=='1' else 'o-Ps','ph':[]}
        elif l.startswith('CREM_REACH'):
            v=l.split()[1:]; Er,Lr=float(v[0]),float(v[1]); e2b,e2a=float(v[5]),float(v[6])
            dem=float(v[7]); h=int(v[10]); nb=float(v[11]); Lb=float(v[12])
            na=nb/math.sqrt(Er); La=Lb*Lr
            cur['ph'].append(dict(nb=nb,Lb=Lb,e2b=e2b,na=na,La=La,e2a=e2a,trim=dem<h-1e-6,dem=dem,h=h))
        elif l.startswith('END'):
            m=re.search(r'stop=(\w+)',l); cur['stop']=m.group(1); cur['fam']=fam.get(cur['i'],'barrier'); recs=[r for r in recs if (r['i'],r['ch'])!=(cur['i'],cur['ch'])]+[cur]
out=open(S+'analysis.txt','w')
def P(*a): print(*a); print(*a,file=out)
P("idx ch   fam      stop    path: n/L/e2 start -> after each photon [T=trimmed by ceiling]")
for r in sorted(recs,key=lambda r:(r['fam'],r['i'],r['ch'])):
    s=" ".join(f"{p['nb']:.3f}/{p['Lb']:.3f}/{p['e2b']:.3f}->{p['na']:.3f}/{p['La']:.3f}/{p['e2a']:.3f}{'T' if p['trim'] else ''}" for p in r['ph'])
    P(f"{r['i']:3d} {r['ch']} {r['fam']:8s} {r['stop']:7s} {s}")
P("\nSUMMARY per family and channel")
for fm in ['circ','ell','barrier']:
  for ch in ['p-Ps','o-Ps']:
    R=[r for r in recs if r['fam']==fm and r['ch']==ch and r['ph']]
    if not R: continue
    last=[r['ph'][-1] for r in R]; first=[r['ph'][0] for r in R]
    P(f"{fm:8s} {ch}: N={len(R)}  photons med {st.median(len(r['ph']) for r in R)}"
      f"  last photon trimmed {sum(p['trim'] for p in last)}/{len(R)}"
      f"  last lands e2<1e-3 {sum(p['e2a']<1e-3 for p in last)}/{len(R)}"
      f"  L_before_last med {st.median(p['Lb'] for p in last):.3f} [{min(p['Lb'] for p in last):.3f},{max(p['Lb'] for p in last):.3f}]"
      f"  L_after_last med {st.median(p['La'] for p in last):.3f}"
      f"  n_after_last med {st.median(p['na'] for p in last):.3f}"
      f"  1-Lb_last vs La: max|diff| {max(abs(abs(1-p['Lb'])-p['La']) for p in last):.2e}")
P("\nCHECK: circular family  n_final = min(L0, 1-L0)  (reflection |L-1| + ceiling landing on the circle)")
for ch in ['p-Ps','o-Ps']:
    R=[r for r in recs if r['fam']=='circ' and r['ch']==ch]
    d=[abs(r['ph'][-1]['na']-min(r['ph'][0]['Lb'],1-r['ph'][0]['Lb'])) for r in R]
    L0=[r['ph'][0]['Lb'] for r in R]
    P(f"  {ch}: N={len(R)} max|n_final - min(L0,1-L0)| = {max(d):.2e}   L0 in [{min(L0):.3f},{max(L0):.3f}]")
P("\nL0 (start, n = 1) per family")
for fm in ['circ','ell','barrier']:
  for ch in ['p-Ps','o-Ps']:
    L0=sorted(r['ph'][0]['Lb'] for r in recs if r['fam']==fm and r['ch']==ch and r['ph'])
    if L0: P(f"  {fm:8s} {ch}: N={len(L0)} L0 = "+" ".join(f"{x:.3f}" for x in L0))
P("\nFIRST photon: harmonic k, target n_k = 1/sqrt(1+2k)?, ceiling circle n = L1 = 1-L0; trimmed iff L1 > n_untrimmed")
for r in sorted(recs,key=lambda r:(r['fam'],r['i'],r['ch'])):
    if not r['ph']: continue
    p=r['ph'][0]
    P(f"  {r['i']:3d} {r['ch']} {r['fam']:8s} k={p['h']} demanded/hbarw={p['dem']:.3f} L1={1-p['Lb']:.3f} n_after={p['na']:.3f} 1/sqrt(1+2k)={1/math.sqrt(1+2*p['h']):.3f} {'T' if p['trim'] else '-'}")

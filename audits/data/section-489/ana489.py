# Audit 489: fractional parts of L/hbar and n_eff right after capture in experiment 5 (CAPDUMP lines).
import math, re, sys
rows=[]
for line in open('exp5-2000.txt'):
    if not line.startswith('CAPDUMP'): continue
    out=re.search(r'out=(\S+)',line).group(1)
    if out not in ('Para-Positronium','Ortho-Positronium'): continue
    L=float(re.search(r' L=(\S+)',line).group(1)); E=float(re.search(r'Ecoul=\S+/(\S+)',line).group(1))
    rows.append((out,L,E))
def ks_uniform(x):
    x=sorted(x); n=len(x); d=max(max((i+1)/n-v, v-i/n) for i,v in enumerate(x)) if n else float('nan')
    # asymptotic p-value
    lam=(math.sqrt(n)+0.12+0.11/math.sqrt(n))*d if n else 0
    p=2*sum((-1)**(k-1)*math.exp(-2*k*k*lam*lam) for k in range(1,100)) if n else float('nan')
    return d,max(0.0,min(1.0,p))
def report(name,sel):
    L=[r[1] for r in sel]; E=[r[2] for r in sel]; n=len(sel)
    if n==0: print(f'{name}: no events'); return
    fL=[v-math.floor(v) for v in L]; neff=[math.sqrt(6.80285/abs(e)) for e in E if e<0]; fn=[v-math.floor(v) for v in neff]
    for lab,f in (('frac(L/hbar)',fL),('frac(n_eff)',fn)):
        d,p=ks_uniform(f); mid=sum(1 for v in f if 0.4<=v<0.6)/len(f); h=[0]*10
        for v in f: h[min(9,int(v*10))]+=1
        print(f'{name} {lab}: n {len(f)}, KS D {d:.3f} p {p:.3f}; fraction in [0.4, 0.6) {mid:.3f} (uniform 0.200 +- {math.sqrt(0.16/len(f)):.3f}); histogram {h}')
    sL=sorted(L); sn=sorted(neff)
    print(f'{name} L/hbar: min {sL[0]:.3f} median {sL[len(sL)//2]:.3f} max {sL[-1]:.3f};  n_eff: min {sn[0]:.3f} median {sn[len(sn)//2]:.3f} max {sn[-1]:.3f}')
print(f'captured events: {len(rows)}')
report('p-Ps',[r for r in rows if r[0].startswith('Para')])
report('o-Ps',[r for r in rows if r[0].startswith('Ortho')])
report('all ',rows)
# H489d: cos(mu_e, mu_p) at capture and averaged over the bound window (CREM_DEBUG_CAPTURE, stderr)
at=[]; mean=[]
for line in open('exp5-2000.err'):
    if line.startswith('CAPTURE '):
        a=re.search(r'atCapture=(\S+)',line).group(1); m=re.search(r'mean=(\S+)',line).group(1)
        try: at.append(float(a)); mean.append(float(m))
        except ValueError: pass
for lab,v in (('cos(mu_e,mu_p) at capture',at),('cos(mu_e,mu_p) bound-window mean',mean)):
    if not v: print(lab,': none'); continue
    u=[(x+1)/2 for x in v]; d,p=ks_uniform(u); n=len(v); h=[0]*10
    for x in v: h[min(9,int((x+1)/0.2))]+=1
    print(f'{lab}: n {n}, KS vs uniform D {d:.3f} p {p:.3f}; histogram (-1..1, 10 bins) {h}')
    for tgt,name in ((1.0,'+-1 (hbar/2)'),(1/math.sqrt(3),'+-0.577 (sqrt3/2)'),(0.5,'+-0.5 (hbar)')):
        w=0.05; f=sum(1 for x in v if abs(abs(x)-tgt)<w)/n
        exp=(min(1,tgt+w)-max(0,tgt-w))   # uniform density 1/2 on [-1,1], both signs -> width total = 2*(band) * 1/2
        print(f'    within 0.05 of {name}: {f:.3f} (uniform {exp:.3f} +- {math.sqrt(exp*(1-exp)/n):.3f})')

# Audit 495: capture in exp. 5 across beam energies -- does capture set l (Langer trace at L = hbar/2)?
# L_in = b sqrt(2 mu K) from each event's own K (centre-of-momentum kinetic energy, eV) and b; L = CAPDUMP L/hbar after the bound window.
# (i) threshold L_c(K) between captures (Ortho/Para) and scatterings, global fit L_c = A (K/1e-4 eV)^-p, per-band thresholds vs the fit;
# (ii) Delta L = L - L_in for captures vs L_in; (iii) f = L mod 1 of captures.  python3 ana495.py <scan dirs...> [--with489 file]
import math, sys, glob, os
hbar=1.054571817e-34; mu=9.1093837015e-31/2; eV=1.602176634e-19
def load(path):
    R=[]
    for l in open(path):
        if not l.startswith('CAPDUMP'): continue
        d=dict(x.split('=',1) for x in l.split()[1:])
        K=float(d['K']); b=float(d['b'])
        e0,e1=(float(x) for x in d['Ecoul'].split('/'))
        out=d['out']
        # energy-gaining "scattering" (final Coulomb-frame energy > 2 x initial; all with rmin < 2 pm): the stiff dipole-barrier region,
        # energy not conserved -- excluded from the threshold as an artefact (audit 495, decided after inspection, reported both ways)
        if out=='Scattering' and e1>2*e0: out='Artefact'
        R.append(dict(out=out,K=K,b=b,L=float(d['L']),Lin=b*math.sqrt(2*mu*K*eV)/hbar,rmin=float(d['rmin'])))
    return R
args=sys.argv[1:]; extra=None
if '--with489' in args: k=args.index('--with489'); extra=args[k+1]; args=args[:k]+args[k+2:]
bands=[]
for d in args:
    R=load(os.path.join(d,'exp5.txt')); bands.append((os.path.basename(d.rstrip('/')),R))
if extra: bands.append(('489',load(extra)))
cap=lambda r: r['out'] in ('Ortho-Positronium','Para-Positronium')
sca=lambda r: r['out']=='Scattering'
allR=[r for _,R in bands for r in R]
print('outcomes per band (capture / scattering / collision / unresolved / artefact), K median, L_in range:')
for name,R in bands:
    c=[r for r in R if cap(r)]; s=[r for r in R if sca(r)]; co=[r for r in R if r['out']=='Collision']; u=[r for r in R if r['out']=='Unresolved']
    Ks=sorted(r['K'] for r in R); art=[r for r in R if r['out']=='Artefact']
    print(f"  {name:10s} n {len(R):4d}: {len(c):3d} / {len(s):4d} / {len(co):3d} / {len(u):3d} / {len(art):3d}"+(f" (rmin max {max(r['rmin'] for r in art)*1e12:.2f} pm)" if art else "")+f"   K_med {Ks[len(Ks)//2]:.3e} eV   L_in {min(r['Lin'] for r in R):.3f}-{max(r['Lin'] for r in R):.3f}")
# L conservation check on scatterings (L_in vs reported L)
sc=[r for r in allR if sca(r)]
dev=max(abs(r['L']-r['Lin'])/r['Lin'] for r in sc)
print(f"check: scatterings |L - L_in|/L_in max {dev:.2e} (n {len(sc)})")
# (i) global separation fit: capture iff L_in (K/1e-4)^p < A
C=[r for r in allR if cap(r)]; S=[r for r in allR if sca(r)]
best=None
for ip in range(0,401):
    p=ip*0.001
    xc=max(math.log(r['Lin'])+p*math.log(r['K']/1e-4) for r in C)
    xs=min(math.log(r['Lin'])+p*math.log(r['K']/1e-4) for r in S)
    mis=sum(1 for r in C if math.log(r['Lin'])+p*math.log(r['K']/1e-4)>xs)+sum(1 for r in S if math.log(r['Lin'])+p*math.log(r['K']/1e-4)<xc)
    gap=xs-xc
    key=(gap>0, gap if gap>0 else -mis)
    if best is None or key>best[0]: best=(key,p,xc,xs)
(_,p,xc,xs)=best
A=math.exp(0.5*(xc+xs))
print(f"(i) global boundary L_c = A (K/1e-4 eV)^-p: p {p:.3f}, A {A:.4f} hbar, log-gap {xs-xc:+.4f} ({'separable' if xs>xc else 'overlap'})")
print("    per band: highest captured L_in, lowest scattered L_in, threshold at band K_med (midpoint), fit at K_med, residual / half-gap")
rows=[]
for name,R in bands:
    c=[r for r in R if cap(r)]; s=[r for r in R if sca(r)]
    if not c or not s: print(f"    {name:10s} -- no bracket ({len(c)} captures, {len(s)} scatterings)"); continue
    Ks=sorted(r['K'] for r in R); Km=Ks[len(Ks)//2]
    # scale each event to the band median K with the global exponent before bracketing
    hc=max(r['Lin']*(r['K']/Km)**p for r in c); ls=min(r['Lin']*(r['K']/Km)**p for r in s)
    mid=0.5*(hc+ls); half=0.5*abs(ls-hc); fit=A*(Km/1e-4)**-p
    rows.append((Km,mid,half,fit))
    print(f"    {name:10s} K_med {Km:.3e}: {hc:.4f} | {ls:.4f}  L_c {mid:.4f} +- {half:.4f}  fit {fit:.4f}  (L_c - fit)/half {((mid-fit)/half if half>0 else float('nan')):+.2f}"
          +("  <-- L_c in [0,45; 0,55]" if 0.45<=mid<=0.55 else ""))
# (ii) Delta L for captures
print("(ii) captures: Delta L = L - L_in vs L_in")
pts=sorted((r['Lin'],r['L']-r['Lin']) for r in C)
for lo,hi in [(0,0.2),(0.2,0.3),(0.3,0.4),(0.4,0.45),(0.45,0.5),(0.5,0.55),(0.55,0.6),(0.6,0.7),(0.7,0.8),(0.8,1.0)]:
    x=[d for L,d in pts if lo<=L<hi]
    if x: print(f"    L_in [{lo:.2f}, {hi:.2f}): n {len(x):3d}  Delta L mean {sum(x)/len(x):+.5f}  min {min(x):+.5f}  max {max(x):+.5f}")
n=len(pts); mx=sum(a for a,_ in pts)/n; my=sum(b for _,b in pts)/n
sxx=sum((a-mx)**2 for a,_ in pts); slope=sum((a-mx)*(b-my) for a,b in pts)/sxx
res=[b-my-slope*(a-mx) for a,b in pts]; se=math.sqrt(sum(e*e for e in res)/(n-2)/sxx)
print(f"    linear fit Delta L = c + s L_in: s {slope:+.5f} +- {se:.5f};  mean Delta L {my:+.5f};  max |Delta L| {max(abs(b) for _,b in pts):.5f}")
# (iii) f of captures
f=sorted(r['L']%1.0 for r in C)
d=max(max((i+1)/n-v, v-i/n) for i,v in enumerate(f)); lam=(math.sqrt(n)+0.12+0.11/math.sqrt(n))*d
pks=max(0.0,min(1.0,2*sum((-1)**(k-1)*math.exp(-2*k*k*lam*lam) for k in range(1,100))))
h=[0]*10
for v in f: h[min(9,int(v*10))]+=1
mid=sum(1 for v in f if 0.4<=v<0.6)/n
print(f"(iii) captures n {n}: L range {min(r['L'] for r in C):.3f}-{max(r['L'] for r in C):.3f}; f histogram (10 bins) {h}; KS vs uniform p {pks:.3g};"
      f" share in [0,4; 0,6) {mid:.3f} (uniform 0,200 +- {math.sqrt(0.2*0.8/n):.3f})")

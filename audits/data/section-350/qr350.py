# Audit 350: annihilation from the n = 1 terminal state with the Quigg-Rosner contact n_QR = (1/2 pi a^3)(hbar/L).
# Input: section-349 out/ files (END lines).  tau_i = 1/Gamma_i, Gamma_i = sigma v [w + (1-w) eps] n_QR.
import glob,os,re,math,statistics as st,sys
D=sys.argv[1]; T0=124.4942e-12  # 1/(sigma v |psi0|^2), engine constants (342)
EPS=1/1113.8978; TP,TO=125.14e-12,142.04e-9
rows=[]
for f in sorted(glob.glob(os.path.join(D,'out','*.txt'))):
    v=os.path.basename(f).split('-')[0]
    for l in open(f):
        if not l.startswith('END'): continue
        x=l.split(); ch='p-Ps' if x[3]=='para' else 'o-Ps'; d=dict(t.split('=') for t in x[4:])
        n=1/float(d['nE']); L=float(d['L']); c=float(d['n/psi0']); G=float(d['Gf'])
        br=G*T0/c if c>0 else float('nan')          # w + (1-w) eps
        w=(br-EPS)/(1-EPS)
        tau=T0*2*L/br if br>0 else float('nan')    # 1/(sigma v br |psi0|^2 hbar/(2L))
        rows.append(dict(v=v,i=int(x[2]),ch=ch,stop=d['stop'],ph=int(d['photons']),n=n,L=L,w=w,br=br,tau=tau))
out=open(os.path.join(D,'qr350.txt'),'w')
def P(*a): print(*a); print(*a,file=out)
P("Quigg-Rosner annihilation from n = 1 (349 terminal states).  T0 = 1/(sigma v |psi0|^2) = 124.4942 ps")
res={}
for v in 'AC':
  for ch in ['p-Ps','o-Ps']:
    for label,keep in [('n=1, no barrier',lambda r:abs(r['n']-1)<1e-6 and r['stop']!='barrier'),
                       ('n=1, all (barrier incl.)',lambda r:abs(r['n']-1)<1e-6)]:
      R=[r for r in rows if r['v']==v and r['ch']==ch and keep(r)]
      if not R: continue
      t=[r['tau'] for r in R]; meas=TP if ch=='p-Ps' else TO
      mean=st.mean(t); harm=1/st.mean(1/x for x in t)
      # survival of the mixture; late slope between 3 and 6 mean lifetimes
      S=lambda tt: st.mean(math.exp(-tt/x) for x in t)
      t1,t2=3*mean,6*mean; late=(t2-t1)/math.log(S(t1)/S(t2))
      P(f"{v}' {ch} [{label}] N={len(R)}: <L>={st.mean(r['L'] for r in R):.4f}  w med {st.median(r['w'] for r in R):.4f} [{min(r['w'] for r in R):.3g},{max(r['w'] for r in R):.3g}]")
      P(f"     mean tau {mean:.4e} s ({mean/meas:.4f} x meas.)  median {st.median(t):.4e} s  1/<Gamma> {harm:.4e} s  <1/G><G> {mean/harm:.3f}  late-slope tau(3-6 mean) {late:.4e} s ({late/meas:.3f} x)  S(mean)={S(mean):.3f}")
      res[(v,ch,label)]=mean
  for label in ['n=1, no barrier','n=1, all (barrier incl.)']:
    if (v,'p-Ps',label) in res and (v,'o-Ps',label) in res:
      P(f"{v}' tau_o/tau_p [{label}] = {res[(v,'o-Ps',label)]/res[(v,'p-Ps',label)]:.1f} (1/eps 1113.9; measured 1135.0)")
P("\nPER TRAJECTORY (n = 1): var idx ch stop L w tau")
for r in sorted(rows,key=lambda r:(r['v'],r['i'],r['ch'])):
    if abs(r['n']-1)<1e-6: P(f"{r['v']} {r['i']:2d} {r['ch']} {r['stop']:7s} L={r['L']:.4f} w={r['w']:.4g} tau={r['tau']:.4e}")

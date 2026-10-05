# Audit 348: cascade from n = 2, variants A/B/C, per channel.
import glob,re,math,statistics as st,os,sys
S=sys.argv[1] if len(sys.argv)>1 else '.'
eV=1.602176634e-19; LINE=5.10179; TP,TO=125.14e-12,142.04e-9; LS=0.17785
recs={}
for f in sorted(glob.glob(os.path.join(S,'out','*.txt'))):
    v,i=os.path.basename(f)[:-4].split('-'); cur=None
    for l in open(f):
        if l.startswith('BEGIN'): _,_,ph=l.split(); cur={'ph':[], 'ch':'p-Ps' if ph=='1' else 'o-Ps'}
        elif l.startswith('CREM_REACH'):
            x=l.split()[1:]; Er=float(x[0]); nb=float(x[11]); Lb=float(x[12])
            cur['ph'].append(dict(E=float(x[4])/eV,k=int(x[10]),trim=float(x[7])<int(x[10])-1e-6,nb=nb,Lb=Lb,
                                  na=nb/math.sqrt(Er),La=Lb*float(x[1]),e2a=float(x[6])))
        elif l.startswith('END'):
            d=dict(t.split('=') for t in l.split()[4:]); cur['stop']=d['stop']
            cur['n']=1/float(d['nE']); cur['L']=float(d['L']); cur['e']=float(d['e'])
            cur['c']=float(d['n/psi0']); cur['ET']=float(d['ET']); recs[(v,int(i),cur['ch'])]=cur
out=open(os.path.join(S,'analysis.txt'),'w')
def P(*a): print(*a); print(*a,file=out)
P("Per variant and channel (n = 2 start, 48 pairs).  Line 2->1 measured 5.10179 eV; L* = 0.17785 hbar (347).")
for v in 'AC':
  for ch in ['p-Ps','o-Ps']:
    R=[r for (vv,i,c),r in sorted(recs.items()) if vv==v and c==ch]
    if not R: continue
    T=TP if ch=='p-Ps' else TO
    first=[r['ph'][0] for r in R if r['ph']]
    ET=[r['ET'] for r in R if math.isfinite(r['ET']) and r['ET']>0]
    stops={}; [stops.__setitem__(r['stop'],stops.get(r['stop'],0)+1) for r in R]
    P(f"\n{v} {ch}: N={len(R)} stops={stops} photons med {st.median(len(r['ph']) for r in R)} max {max(len(r['ph']) for r in R)}")
    if first:
      P(f"   first photon: E/eV med {st.median(p['E'] for p in first):.4f}; at 5.10179+-0.01: {sum(abs(p['E']-LINE)<0.01 for p in first)}/{len(first)};"
        f" k=3 untrimmed {sum(p['k']==3 and not p['trim'] for p in first)}; trimmed {sum(p['trim'] for p in first)}/{len(first)}")
      hist={}; [hist.__setitem__(round(p['E'],2),hist.get(round(p['E'],2),0)+1) for p in first if not p['trim']]
      P("   untrimmed first-photon energies (eV:count): "+" ".join(f"{k}:{c}" for k,c in sorted(hist.items())))
    P(f"   terminal n: med {st.median(r['n'] for r in R):.4f} [{min(r['n'] for r in R):.3f},{max(r['n'] for r in R):.3f}];"
      f" n in [0.95,1.05]: {sum(0.95<=r['n']<=1.05 for r in R)}/{len(R)}")
    P(f"   terminal L: med {st.median(r['L'] for r in R):.4f} [{min(r['L'] for r in R):.3f},{max(r['L'] for r in R):.3f}];"
      f" in L*+-10% [0.160,0.196]: {sum(0.160<=r['L']<=0.196 for r in R)}/{len(R)}; e med {st.median(r['e'] for r in R):.3f}")
    P(f"   contact n/psi0: med {st.median(r['c'] for r in R):.3e} [{min(r['c'] for r in R):.2e},{max(r['c'] for r in R):.2e}]")
    if ET: P(f"   E[T]: N={len(ET)} mean {st.mean(ET):.4e} s ({st.mean(ET)/T:.3g} x meas.), median {st.median(ET):.4e} s ({st.median(ET)/T:.3g} x)")
  same=sum(1 for (vv,i,c),r in recs.items() if vv==v and c=='p-Ps' and (v,i,'o-Ps') in recs
           and recs[(v,i,'o-Ps')]['stop']==r['stop'] and len(recs[(v,i,'o-Ps')]['ph'])==len(r['ph']))
  P(f"   {v}: pairs with same stop & photon count: {same}/{sum(1 for k in recs if k[0]==v and k[2]=='p-Ps')}")
  mp=[recs[(v,i,'p-Ps')]['ET'] for (vv,i,c) in recs if vv==v and c=='p-Ps' and (v,i,'o-Ps') in recs]
  mo=[recs[(v,i,'o-Ps')]['ET'] for (vv,i,c) in recs if vv==v and c=='p-Ps' and (v,i,'o-Ps') in recs]
  ok=[(a,b) for a,b in zip(mp,mo) if math.isfinite(a) and math.isfinite(b) and a>0 and b>0]
  if ok: P(f"   {v}: tau_o/tau_p (means) = {st.mean(b for a,b in ok)/st.mean(a for a,b in ok):.1f} (measured 1135)")
P("\nACTION-RULE CHECKS (349)")
for v in 'AC':
  for ch in ['p-Ps','o-Ps']:
    R=[r for (vv,i,c),r in sorted(recs.items()) if vv==v and c==ch]
    if not R: continue
    E=[p['E'] for r in R for p in r['ph']]
    P(f"{v} {ch}: photons {len(E)}; E/eV min {min(E) if E else float('nan'):.6f} max {max(E) if E else float('nan'):.6f}; vs 5.10212: max|d| {max(abs(x-5.10212) for x in E) if E else float('nan'):.2e}; vs measured 5.10179: {(st.mean(E)-LINE)*1e3 if E else float('nan'):+.3f} meV")
    nb=[r for r in R if r['stop']!='barrier' and r['ph']]
    P(f"   non-barrier with photons: {len(nb)}; max|n-1| {max(abs(r['n']-1) for r in nb) if nb else float('nan'):.2e};"
      f" max|L - |L0-1|| {max(abs(r['L']-abs(r['ph'][0]['Lb']-1)) for r in nb) if nb else float('nan'):.2e};"
      f" max|J_r change| (n-L before vs after) {max(abs((r['n']-r['L'])-(r['ph'][0]['nb']-r['ph'][0]['Lb'])) for r in nb) if nb else float('nan'):.2e}")
    z=[r for r in R if not r['ph']]
    P(f"   zero-photon runs: {len(z)}; stops {sorted(set(r['stop'] for r in z))}; n {sorted(set(round(r['n'],4) for r in z))}")
    L0=[r['ph'][0]['Lb'] for r in R if r['ph']]
    P(f"   L0 of emitting runs: min {min(L0) if L0 else float('nan'):.3f}; L0<1 among emitting: {sum(x<1 for x in L0)}")
    bar=[r for r in R if r['stop']=='barrier']
    P(f"   barrier: {len(bar)}, terminal L {[round(r['L'],4) for r in bar]}")
P("\nPER TRAJECTORY: var idx ch stop | photons E/eV(k,T) | terminal n L e contact ET")
for (v,i,c),r in sorted(recs.items()):
    P(f"{v} {i:2d} {c} {r['stop']:7s} | "+" ".join(f"{p['E']:.3f}({p['k']}{'T' if p['trim'] else ''})" for p in r['ph'])
      +f" | n={r['n']:.4f} L={r['L']:.4f} e={r['e']:.4f} c={r['c']:.3e} ET={r['ET']:.3e}")

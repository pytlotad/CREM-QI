import re,statistics as st,math
S='/tmp/claude-1000/-home-teddy-Projekty-Positronium/5981d9ee-dfbe-485a-8f93-5c6d87fd7653/scratchpad/sec343/'
rows=[]
for l in open(S+'final.txt'):
    m=re.match(r'idx\s+(\d+) (\w+)\s+stop=(\w+) photons=(\d+) (.*)',l)
    if not m: continue
    d={k:float(v) for k,v in (x.split('=') for x in m.group(5).split())}
    d['nE']=1.0/d['nE']   # the probe printed sqrt(a1/A) = 1/n_E
    d.update(idx=int(m.group(1)),ch=m.group(2),stop=m.group(3),ph=int(m.group(4)))
    if math.isfinite(d['ET']): rows.append(d)
fmt=lambda v:f"{v:.3g}"
for ch in ('para','ortho'):
    R=[r for r in rows if r['ch']==ch]
    fam={'kolowe (e<0,05)':[r for r in R if r['stop']=='closed' and r['e']<0.05],
         'eliptyczne':[r for r in R if r['stop']=='closed' and r['e']>=0.05],
         'bariera':[r for r in R if r['stop']=='barrier']}
    ET=[r['ET'] for r in R]
    print(f"== {ch}: n={len(R)}, mean E[T] {st.mean(ET)*1e12:.4g} ps")
    for name,F in fam.items():
        if not F: print(f"   {name:16s}: 0"); continue
        et=[r['ET'] for r in F]; share=sum(et)/sum(ET)
        print(f"   {name:16s}: {len(F):2d}  nE med {st.median(r['nE'] for r in F):.3f} [{min(r['nE'] for r in F):.3f},{max(r['nE'] for r in F):.3f}]"
              f"  L med {st.median(r['L'] for r in F):.3f}  e med {st.median(r['e'] for r in F):.3f}"
              f"  rp/eps med {st.median(r['rp/eps'] for r in F):.3g}  n/psi0 med {st.median(r['n/psi0'] for r in F):.3g}"
              f"  E[T] med {st.median(et)*1e12:.4g} ps, share of ensemble mean {100*share:.1f}%")
    C=[r for r in R if r['stop']=='closed']
    print(f"   closure margin |L-1| - nE: min {min(abs(r['L']-1)-r['nE'] for r in C):+.3e}")

# Audit 478: spin rotation in one close pass vs r_min (mechanical path).
import math, sys
ONLYESC=int(sys.argv[1]) if len(sys.argv)>1 else 0
rows=[]
for l in open('scan.txt'):
    t=l.split()
    if t[0]!='b': continue
    rows.append(dict(b=float(t[1]),rmin=float(t[4]),beta=float(t[7]),out=int(t[9]),stop=int(t[11]),c0=float(t[13]),cf=float(t[15]),
                     rot1=float(t[19]),rot2=float(t[21]),sep=float(t[24])))
re=2.8179403262
print('b[fm]  n  <r_min>   beta_max  escaped  <rot1>    <rot2>    <rot> r_min/r_e  <|dcos12|>/<rot>  frac rot>2.5  histogram rot (0..pi, 8 bins)')
for b in sorted(set(r['b'] for r in rows),reverse=True):
    rs=[r for r in rows if r['b']==b and (ONLYESC==0 or r['sep']>3000)]; n=len(rs)
    if n==0: continue
    rm=sum(r['rmin'] for r in rs)/n; be=max(r['beta'] for r in rs); esc=sum(1 for r in rs if r['sep']>3000)
    rot=[r['rot1'] for r in rs]+[r['rot2'] for r in rs]; mr=sum(rot)/len(rot)
    dc=sum(abs(r['cf']-r['c0']) for r in rs)/n
    big=sum(1 for x in rot if x>2.5)/len(rot)
    h=[0]*8
    for x in rot: h[min(7,int(x/math.pi*8))]+=1
    print(f'{b:6.0f} {n:3d} {rm:9.3f} {be:8.4f} {esc:4d}/{n}  {sum(r["rot1"] for r in rs)/n:.3e} {sum(r["rot2"] for r in rs)/n:.3e}  {mr*rm/re:8.3f}   {dc/mr if mr>0 else 0:8.3f}       {big:.3f}     {h}')

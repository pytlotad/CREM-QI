# Audit 613 T4: offline G shares (var599 rate, both normalized within Delta l = -1) vs engine p613.txt
import re
s599=open('../section-599/var599.py').read().split('rows=[]')[0]; g={}; exec(s599,g)
worst=0
for line in open('p613.txt'):
    head,rest=line.split(':'); n=int(head[0]); l='spdf'.index(head[1])
    eng={int(a):float(b) for a,b in re.findall(r"n'=(\d+) ([0-9.e-]+)",rest)}
    r={n2:g['rate']('G','o',n,l,n2,l-1) for n2 in eng}; t=sum(r.values())
    d=max(abs(eng[k]-r[k]/t) for k in eng); worst=max(worst,d); print(head, f"max |engine - offline| = {d:.1e}")
print(f"T4 max {worst:.1e}")

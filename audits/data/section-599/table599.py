# Audit 599 (report after the verdict): lifetimes tau/QM state by state, one place, variants Ao (today's midpoint) and Go (geometric).
exec(open('var599.py').read().split('rows=[]')[0])
print("state   Ao (one place)   Go (one place)")
for l in range(4):
    for n in range(max(l+1,2),9):
        if (n,l)==(2,0): continue
        q=sum(rate_qm(n,l,n2,l2) for n2,l2 in channels(n,l))
        out=[]
        for k in ('A','G'):
            r=sum(rate(k,'o',n,l,n2,l2) for n2,l2 in channels(n,l)); out.append(q/r)
        print(f"{n}{'spdf'[l]}      {out[0]:6.3f}          {out[1]:6.3f}")

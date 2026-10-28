import re,math,sys
for l in open(sys.argv[1]):
    if not l.startswith("CREM_SPINVEC"): continue
    d=dict(re.findall(r"(\w+)=([^ ]+)",l)); v=lambda k:[float(x) for x in d[k].split(",")]
    m2,L=v("m2"),v("L"); n=lambda a: math.sqrt(sum(x*x for x in a))
    print("%.3e %+.6f"%(float(d["t"]),-sum(x*y for x,y in zip(L,m2))/n(L)/n(m2)),end="; ")
print()

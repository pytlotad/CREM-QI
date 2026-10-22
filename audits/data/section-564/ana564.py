#!/usr/bin/env python3
# Audit 564: engine F1+F2 against the field-momentum expectation, per separation, speed and R.
import glob, os, re, sys
d=sys.argv[1] if len(sys.argv)>1 else '.'
rows=[]
for f in sorted(glob.glob(os.path.join(d,'sep1','b*_R*.txt'))+glob.glob(os.path.join(d,'sepN','s*_R*.txt'))):
    t=open(f).read()
    m=re.search(r'beta (\S+)\s+nodes (\d+)\s+span \S+ r/c\s+R (\S+) r\*',t)
    sep=re.search(r'separation (\S+) r\*',t); sep=float(sep.group(1)) if sep else 1.0
    net=[float(x) for x in re.search(r'engine F1\+F2\s+(\S+) (\S+) (\S+)',t).groups()]
    exp=[float(x) for x in re.search(r'expected\s+(\S+) (\S+) (\S+)',t).groups()]
    F1=re.search(r'\(\|F1\| (\S+) N',t)
    st=float(re.search(r'\|P_V \+ p_hid\|/\|p_hid\|\s+(\S+)',t).group(1))
    rows.append((sep,float(m.group(1)),float(m.group(3))/sep,int(m.group(2)),net,exp,st,os.path.basename(f)))
print(f"{'sep':>5} {'beta':>5} {'R/sep':>5} {'nodes':>5} {'|F1+F2| N':>11} {'|expected| N':>12} {'exp.x/net.x':>11} {'|net-exp|/|net|':>15} {'|P_V+p_h|/|p_h|':>15}")
for sep,b,R,n,net,exp,st,f in sorted(rows):
    nn=sum(x*x for x in net)**.5; ne=sum(x*x for x in exp)**.5
    dd=sum((a-c)**2 for a,c in zip(net,exp))**.5
    print(f"{sep:5.0f} {b:5.2f} {R:5.1f} {n:5d} {nn:11.4e} {ne:12.4e} {exp[0]/net[0]:11.4f} {dd/nn:15.4f} {st:15.4e}")

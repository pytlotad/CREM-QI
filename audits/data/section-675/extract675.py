# Audit 675: per-photon table from 668 probes for img675.cpp.  usage: python3 extract675.py ../section-668/runs > photons675.txt
import re,glob,sys
num=r'(\S+)'
def v3(t,key): m=re.search(key+' '+' '.join([num]*3),t); return m.groups()
for fn in sorted(glob.glob(sys.argv[1]+'/A*_*.err')):
    t=open(fn).read()
    d=re.search(r'PAIR eigen draw .* a1 (\S+) a2 (\S+) K (\S+) w_k \S+ E1S (\S+)',t)
    lg=re.search(r'PAIR eigen ledger .* after (\S+) B \S+ S2P \S+ S1S (\S+) Zsb \S+ Zsa (\S+)',t)
    q=re.search(r'PAIR field eigen \(1S\) q=\S+ m=(-?\d+) P_S=(\S+)',t)
    if not (d and lg and q): continue
    Bs=t.split('| B ')[1].split('| A ')[0]; As=t.split('| A ')[1]
    EA=re.search(r'^E (\S+) L (\S+)',As).groups()
    sf=re.search(r'CREM_LEDGER t \S+ Eg \S+ harm \S+ hel \S+ singlet (\d)',t).group(1)
    row=[fn.split('/')[-1],sf]+list(v3(t,'BHAT'))+list(d.groups())+list(lg.groups())+list(q.groups())+list(EA)+list(v3(As,'ax'))+list(v3(Bs,'peri'))+list(v3(As,'m1'))+list(v3(As,'m2'))
    print(' '.join(row))

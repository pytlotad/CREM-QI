# Audit 517: S band -- new (k = 1) vs import 17 vs bare E_rel; check (1) ladder, (2) Delta E = B (1 - 1/(2L)), B = E_rel(1/2) - E_rel(1)
# (bare), (3) jump at L = hbar, 2P1/2 (import 21, lambda = g - 1 and lambda = 1).  python3 ana517.py
def load(f):
    F={}; P={}
    for l in open(f):
        t=l.split()
        if t[0]=='F': F[(t[1],int(t[2]),float(t[3]))]=float(t[4])
        else: P[(t[1],int(t[2]))]=float(t[3])
    return F,P
N,PN=load('p517-default.txt'); O,_=load('p517-import17.txt'); R,_=load('p517-bare.txt')
pairs=['proton,electron','antimuon,electron','electron,positron']
for p in pairs:
  for n in (1,2,3,5):
    B=R[(p,n,0.5)]-R[(p,n,1.0)]
    ks=sorted(L for (q,m,L) in N if q==p and m==n)
    lad=N[(p,n,0.5)]-O[(p,n,0.5)]
    dev=max(abs((N[(p,n,L)]-O[(p,n,L)])-B*(1-1/(2*L))) for L in ks if L<1)
    row=' '.join(f"{L:g}:{N[(p,n,L)]-O[(p,n,L)]:+.1f}" for L in ks if L<1 and L!=0.999999)
    jo=O[(p,n,0.999999)]-O[(p,n,1.000001)] if (p,n,1.000001) in O else float('nan')
    jn=N[(p,n,0.999999)]-N[(p,n,1.000001)] if (p,n,1.000001) in N else float('nan')
    print(f"{p:18s} n {n}: B {B:+.3f} MHz; ladder new - old {lad:+.2e} MHz; new - old (MHz) {row}; max |new - old - B(1 - 1/2L)| {dev:.2e};"
          f" F(hbar-) - F(hbar+): import 17 {jo:+.3f} (-B/2 {-B/2:+.3f}), new {jn:+.2e}")
    if (p,n) in PN:
      # 2P1/2 with lambda = g - 1 (default) ; with lambda = 1 it is F(k=1) = F(hbar)
      print(f"{'':18s}       S-band edge F(hbar-) new {N[(p,n,0.999999)]:+.4f} = 2P1/2 at lambda = 1 F(hbar) {R[(p,n,1.0)]:+.4f}; 2P1/2 engine (lambda = g - 1) {PN[(p,n)]:+.4f}"
            f" -- edge - 2P1/2 {N[(p,n,0.999999)]-PN[(p,n)]:+.4f} MHz (the anomalous part of 2P1/2); import 17 edge - 2P1/2 {O[(p,n,0.999999)]-PN[(p,n)]:+.3f}")

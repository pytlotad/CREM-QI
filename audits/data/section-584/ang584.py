# Audit 584: the weight of the midpoint mode rule (560) against QM -- angular factor l_>/(2l+1), sum rule, all l.
import math, sys, io, contextlib
src=open('../section-560/tempo_e1_punkt_srodkowy.py').read().split('if __name__ == "__main__":')[0]
g={}; exec(src,g)
rate_qm, rate_rule = g['rate_qm'], g['rate_rule']
names="spdfghik"
print("T1  ns -> n'p: A(RB, w = 1)/QM")
for n in range(3,11):
    row=[f"{rate_rule(n,0,n2,1,'RB')/rate_qm(n,0,n2,1):.3f}" for n2 in range(2,n)]
    print(f"  {n}s -> n'p (n' = 2..{n-1}): "+' '.join(row))
print("\nT2  Delta l = +1 from l >= 1 (n' = n - 1 first): A(RBw)/QM  [A(RB)/QM]  w = 2 l_>/(2l+1)")
better=0; total=0
for l in range(1,5):
    for n in range(l+3,9):
        cells=[]
        for n2 in range(n-1,l+1,-1):   # l' = l+1 needs n2 >= l+2
            q=rate_qm(n,l,n2,l+1); rw=rate_rule(n,l,n2,l+1,'RBw'); rb=rate_rule(n,l,n2,l+1,'RB')
            cells.append(f"{n2}{names[l+1]} {rw/q:.3f} [{rb/q:.3f}]")
            if n2==n-1:
                total+=1; better+= abs(math.log(rw/q))<abs(math.log(rb/q))
        print(f"  {n}{names[l]} -> "+'  '.join(cells))
print(f"  n' = n-1 channels where the weight helps: {better}/{total}")
print("\nT3  Delta l = -1 at the midpoint, w = 2l/(2l+1): A(RBw)/QM by Delta n")
for l in range(1,4):
    for n in range(l+2,9):
        cells=[f"{n2}{names[l-1]} {rate_rule(n,l,n2,l-1,'RBw')/rate_qm(n,l,n2,l-1):.3f}" for n2 in range(n-1,l-1,-1)]
        print(f"  {n}{names[l]} -> "+'  '.join(cells))
print("\nsum rule of the QM angular factor: sum_{l'} l_>/(2l+1) = "+', '.join(f"l={l}: {((l+1)+(l if l>0 else 0))/(2*l+1):.3f}" for l in range(0,5)))

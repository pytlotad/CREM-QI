# Audit 585: offline up-share (variant H) for the engine test states, more digits.
src=open('hyb585.py').read().split('ns=lambda')[0]
g={}; exec(src,g)
for n,l in [(4,1),(5,2),(6,1)]:
    AE=g['e_all'](n,l); ch=list(g['channels'](n,l))
    ru=sum(g['rate_rule'](n,l,n2,l2,'RBw') for n2,l2 in ch if l2==l+1)
    q=sum(g['rate_qm'](n,l,n2,l2) for n2,l2 in ch)
    print(f"{n}{'spdf'[l]}  up-share H {ru/(AE+ru):.6f}  tau_H/QM {q/(AE+ru):.4f}  tau_E/QM {q/AE:.4f}")

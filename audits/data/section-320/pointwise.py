exec(open(__file__.replace('pointwise.py','plane_qed.py')).read().split('rng=np.random')[0])
rng=np.random.default_rng(2); rows=[]
while len(rows)<300:
    x1,x2=rng.uniform(0,1,2); x3=2-x1-x2
    if not (0<x3<1): continue
    c12=(x3**2-x1**2-x2**2)/(2*x1*x2)
    k1=x1*np.array([1.,0,0]); k2=x2*np.array([c12,np.sqrt(max(0,1-c12**2)),0]); ks=[k1,k2,-(k1+k2)]
    T=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            e=np.zeros(3); e[i]=1; f=np.zeros(3); f[j]=1
            T[i,j]=0.5*(rate(ks,e+f)-rate(ks,e)-rate(ks,f)) if i!=j else rate(ks,e)
    ev=np.linalg.eigvalsh(T[:2,:2]); rows.append((T[2,2]/(0.5*np.trace(T[:2,:2])), (ev[1]-ev[0])/ev.sum(), abs(T[0,2])+abs(T[1,2])))
r=np.array(rows)
print("pointwise b/a: min %.4f max %.4f"%(r[:,0].min(),r[:,0].max()))
print("in-plane anisotropy (l1-l2)/(l1+l2): min %.4f max %.4f"%(r[:,1].min(),r[:,1].max()))
print("max |T_xz|+|T_yz|: %.2e"%r[:,2].max())

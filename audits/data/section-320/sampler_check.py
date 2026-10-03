import numpy as np
rng=np.random.default_rng(5); N=4_000_000
def f(x):
    u=1-x; d=2-x
    return x*u/d**2-2*u*u*np.log(u)/d**3+d/x+2*u*np.log(u)/x**2
# exact joint Ore-Powell
x1=rng.uniform(0,1,N); x2=rng.uniform(0,1,N); x3=2-x1-x2; ok=(x3>0)&(x3<1)
x1,x2,x3=x1[ok],x2[ok],x3[ok]
M=((1-x1)/(x2*x3))**2+((1-x2)/(x1*x3))**2+((1-x3)/(x1*x2))**2
print("max |M|^2 on region %.4f"%M.max())
acc=rng.uniform(0,2,len(M))<M
e=x1[acc]; print("exact joint:   <x1> %.4f  <x1^2> %.4f  P(x1>0.9) %.4f  <min x> %.4f"%(e.mean(),(e**2).mean(),(e>0.9).mean(),np.minimum(np.minimum(x1,x2),x3)[acc].mean()))
# existing sampler: uniform in region, accept on f(x1)
xx=np.linspace(1e-4,0.9999,10000); b=max(f(0.999),1.0)
a2=rng.uniform(0,1,len(x1))*b<=f(x1)
e=x1[a2]; print("code sampler:  <x1> %.4f  <x1^2> %.4f  P(x1>0.9) %.4f  <min x> %.4f"%(e.mean(),(e**2).mean(),(e>0.9).mean(),np.minimum(np.minimum(x1,x2),x3)[a2].mean()))
# marginal of f alone
w=f(xx); w/=w.sum(); print("spectrum f:    <x>  %.4f  <x^2>  %.4f  P(x>0.9)  %.4f"%((xx*w).sum(),(xx*xx*w).sum(),w[xx>0.9].sum()))

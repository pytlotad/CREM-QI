// Audit 432: hydrogen E1 rates -- exact QM (numerical radial integrals) vs classical mode correspondence variants.
#include <cmath>
#include <cstdio>
#include <vector>
#include <functional>
const double alpha=7.2973525693e-3, tau0=2.4188843265857e-17;   // a.u. of time [s]
double fact(int n){ double f=1; for(int i=2;i<=n;++i) f*=i; return f; }
// hydrogen radial R_nl(r), a.u.
double R(int n,int l,double r){
  const double rho=2*r/n; const int k=n-l-1; const double a=2*l+1;
  // generalized Laguerre L_k^{a}(rho) by recursion
  double L0=1, L1=1+a-rho, Lk=(k==0?L0:L1);
  for(int i=1;i<k;++i){ const double L2=((2*i+1+a-rho)*L1-(i+a)*L0)/(i+1); L0=L1; L1=L2; Lk=L2; }
  const double norm=std::sqrt(std::pow(2.0/n,3)*fact(n-l-1)/(2.0*n*fact(n+l)));
  return norm*std::exp(-rho/2)*std::pow(rho,l)*Lk;
}
double radial(int n,int l,int n2,int l2){
  const double rmax=12.0*std::max(n,n2)*std::max(n,n2)+40; const int N=200000; const double h=rmax/N; double s=0;
  for(int i=0;i<=N;++i){ const double r=i*h, w=(i==0||i==N)?1:(i%2?4:2); s+=w*R(n,l,r)*R(n2,l2,r)*r*r*r; }
  return s*h/3; }
double Aqm(int n,int l,int n2,int l2){ const double w=0.5/(n2*n2)-0.5/(n*n), I=radial(n,l,n2,l2);
  return (4.0/3.0)*std::pow(alpha,3)*w*w*w*(std::max(l,l2)/(2.0*l+1))*I*I/tau0; }
// J_{k-1}, J_k, J_{k+1} at x (x < ~k) by Miller backward recurrence, normalized by J0 + 2 sum J_{2i} = 1
void bessel3(int k,double x,double& jm,double& j0,double& jp){
  const int N=k+20+int(10*std::sqrt(double(k)))+(int)x; std::vector<double> J(N+2,0.0); J[N+1]=0; J[N]=1e-300;
  for(int i=N;i>=1;--i){ J[i-1]=2.0*i/x*J[i]-J[i+1]; if(std::abs(J[i-1])>1e250){ for(int t=i-1;t<=N+1;++t) J[t]*=1e-250; } }
  double norm=J[0]; for(int i=2;i<=N;i+=2) norm+=2*J[i];
  jm=J[k-1]/norm; j0=J[k]/norm; jp=J[k+1]/norm; }
double Jp(int k,double x){ double a,b,c; bessel3(k,x,a,b,c); return 0.5*(a-c); }
double Jk(int k,double x){ double a,b,c; bessel3(k,x,a,b,c); return b; }
// classical mode amplitude |c_{k,m}| (units of a), m = +1 co-rotating, -1 counter
double cmode(int k,double e,int m){ const double x=k*e; const double A=(2.0/k)*Jp(k,x), B=2*std::sqrt(1-e*e)*Jk(k,x)/x; return 0.5*std::abs(A+m*B); }
// classical rate for transition, variant v
double Acl(int v,int n,int l,int n2,int l2){
  const int k=n-n2, m=(l2==l-1)?+1:-1, lg=std::max(l,l2);
  double nb, Lb; if(v==0){ nb=n; Lb=l+0.5; } else if(v==1){ nb=0.5*(n+n2); Lb=lg; } else { nb=n; Lb=lg; }
  const double e=std::sqrt(std::max(0.0,1-Lb*Lb/(nb*nb))), a=nb*nb, c=cmode(k,e,m);
  if(v==0){ const double Om=1.0/(n*n*n); const double P=(2.0/3.0)*std::pow(alpha,3)*std::pow(k*Om,4)*a*a*c*c; return P/(k*Om)/tau0; }
  const double w=0.5/(n2*n2)-0.5/(n*n);
  return (2.0/3.0)*std::pow(alpha,3)*w*w*w*a*a*c*c/tau0; }
// engine-style photon count at the initial Langer orbit: Sum_k P_{k,m}/(k hbar Omega), k <= kmax (kmax < 0: all)
double Aeng(int n,int l,int m,int kmax){
  const double Lb=l+0.5, e=std::sqrt(std::max(0.0,1-Lb*Lb/(double(n)*n))), a=double(n)*n, Om=1.0/(double(n)*n*n);
  auto term=[&](int k){ const double c=cmode(k,e,m); return (2.0/3.0)*std::pow(alpha,3)*std::pow(k*Om,4)*a*a*c*c/(k*Om); };
  double s=0;
  if(kmax>0){ for(int k=1;k<=kmax;++k) s+=term(k); return s/tau0; }
  // all harmonics: exact up to k = 64, then sampled every dk = k/64 (trapezoid), stop when the tail is negligible
  int k=1; for(;k<=64;++k) s+=term(k);
  double prev=term(64); int kp=64;
  while(true){ const int dk=std::max(1,kp/64); const int kn=kp+dk; const double t=term(kn); s+=0.5*(prev+t)*dk-0.0; 
    // trapezoid counts the left end twice at the first step: correct by removing half of term(64) once
    if(kp==64) s-=0.5*prev;
    prev=t; kp=kn; if(t*kn<1e-9*s||kn>2000000) break; }
  return s/tau0; }
int main(){
  std::printf("check A(2p->1s) = %.4e /s (6.2649e8)\n",Aqm(2,1,1,0));
  const int st[13][2]={{3,0},{4,0},{5,0},{6,0},{2,1},{3,1},{4,1},{5,1},{3,2},{4,2},{5,2},{4,3},{5,3}};
  const char* L="spdf";
  for(auto& s: st){ const int n=s[0], l=s[1];
    double tq=0, tc[3]={0,0,0}; std::vector<double> bq, bc;
    for(int n2=1;n2<n;++n2) for(int l2: {l-1,l+1}){ if(l2<0||l2>=n2) continue;
      const double a=Aqm(n,l,n2,l2); tq+=a; bq.push_back(a);
      for(int v=0;v<3;++v){ const double c=Acl(v,n,l,n2,l2); tc[v]+=c; if(v==1) bc.push_back(c); } }
    std::printf("%d%c: tau QM %8.3f ns | K0 %8.3f (x%.3f) | K1 %8.3f (x%.3f) | K3 %8.3f (x%.3f) | branches QM/K1:",n,L[l],1e9/tq,1e9/tc[0],tq/tc[0],1e9/tc[1],tq/tc[1],1e9/tc[2],tq/tc[2]);
    int i=0; for(int n2=1;n2<n;++n2) for(int l2: {l-1,l+1}){ if(l2<0||l2>=n2) continue; std::printf(" %d%c %.3f/%.3f",n2,L[l2],bq[i]/tq,bc[i]/tc[1]); ++i; }
    std::printf("\n");
    // engine-style variants
    if(l==0){ const double net=Aeng(n,0,+1,-1)-Aeng(n,0,-1,-1);
      std::printf("    E_net (engine today, Sum (P+ - P-)/(k hbar Omega)) %9.3f ns (x%.3f)\n",1e9/net,tq/net);
      const double allm=Aeng(n,0,-1,-1), allowm=Aeng(n,0,-1,n-2), allb=allm+Aeng(n,0,+1,-1), allowb=allowm+Aeng(n,0,+1,n-2);
      std::printf("    E_all(m=-1) %9.3f ns (x%.3f)  E_allow(m=-1) %9.3f (x%.3f)  E_all(both) %9.3f (x%.3f)  E_allow(both) %9.3f (x%.3f)\n",
        1e9/allm,tq/allm,1e9/allowm,tq/allowm,1e9/allb,tq/allb,1e9/allowb,tq/allowb); }
    else { const double allp=Aeng(n,l,+1,-1), allowp=Aeng(n,l,+1,n-l);
      std::printf("    E_all(m=+1) %9.3f ns (x%.3f)  E_allow(m=+1, k<=n-l) %9.3f (x%.3f)\n",1e9/allp,tq/allp,1e9/allowp,tq/allowp); }
  }
}

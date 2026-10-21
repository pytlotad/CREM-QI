// Audit 543 (from m542.cpp; basin 2 = on the conserved axis J^, basin 3 = J^ for x_J < 1, s^ for x_J >= 1):
// Audit 542: two-spin + orbit classical model of the Ps 2^3P transport in a field (variants C, A, U1, U2 of prereg542), isotropic 3P_J ensembles,
// singlet fraction at the 2P decay P_S = Gamma int e^{-Gamma t} w dt, w = (1 - n1.n2)/2.  Units: hbar = 1, angular frequencies in rad/s.
#include <cstdio>
#include <cmath>
#include <random>
#include <array>
#include <vector>
#include <string>
using V=std::array<double,3>;
static V operator+(V a,V b){return {a[0]+b[0],a[1]+b[1],a[2]+b[2]};} static V operator-(V a,V b){return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
static V operator*(V a,double s){return {a[0]*s,a[1]*s,a[2]*s};} static double dot(V a,V b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
static V cross(V a,V b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};} static double nrm(V a){return std::sqrt(dot(a,a));}
const double PI=3.14159265358979323846, h=6.62607015e-34, hbar=h/(2*PI), muB=9.2740100783e-24, g=2.0023193044, tau=3.348e-9, G=1/tau;
const double EJ[3]={-7.312222e9*2*PI,-1.826996e9*2*PI,2.558642e9*2*PI}, XJ[3]={0.5,1.5,2.5};   // rad/s; E_S = 0
static void parab(double x,double& H,double& dH){ const double d01=(EJ[1]-EJ[0])/(XJ[1]-XJ[0]), d12=(EJ[2]-EJ[1])/(XJ[2]-XJ[1]), cu=(d12-d01)/(XJ[2]-XJ[0]);
  H=EJ[0]+d01*(x-XJ[0])+cu*(x-XJ[0])*(x-XJ[1]); dH=d01+cu*((x-XJ[0])+(x-XJ[1])); }
struct St{V n1,n2,l;};
static int variant; static double zeta; static int basin=0;   // zeta = Z/hbar
static void deriv(const St& s,St& d){
  const V Je=s.l+(s.n1+s.n2)*0.5; const double x=nrm(Je); V Jh=Je*(1/x); double H,dH; parab(x,H,dH);
  const double c12=dot(s.n1,s.n2); V o1{},o2{},ol{};
  if(variant==0||variant==1||(variant==2&&c12<0)){ if(c12>=0){ o1=Jh*dH; o2=o1; ol=o1; } }   // 530 / 534
  else if(variant==4||variant==5){ const V sm=s.n1+s.n2; const double sn=nrm(sm); const V sh=sm*(1/sn); const V Jp=s.l+sh; const double xp=nrm(Jp); const V Jph=Jp*(1/xp);
    double Hp,dHp; parab(xp,Hp,dHp); const double wT=(3+c12)/4; const V perp=(Jph-sh*dot(Jph,sh))*(2/sn);   // 2/sn: s = sn/2
    o1=perp*(wT*dHp)+s.n2*(0.5*Hp); o2=perp*(wT*dHp)+s.n1*(0.5*Hp); ol=Jph*(wT*dHp); }
  else { const double wT=(3+c12)/4; o1=Jh*(wT*dH)+s.n2*(0.5*H); o2=Jh*(wT*dH)+s.n1*(0.5*H); ol=Jh*(wT*dH); }
  const V z{0,0,1};
  if(variant==0||variant==5){ o1=o1+z*(2*zeta); o2=o2-z*(2*zeta); }
  else { const V sum=s.n1+s.n2; const double sm=0.5*nrm(sum);
    if(sm>1e-6){ const double U=zeta*(s.n1[2]-s.n2[2]); const V cm=sum*(U/(4*sm)); o1=o1+(z*(sm*zeta)+cm)*2; o2=o2+(z*(-sm*zeta)+cm)*2; }
    else { o1=o1+z*(2*zeta); o2=o2-z*(2*zeta); } }
  d.n1=cross(o1,s.n1); d.n2=cross(o2,s.n2); d.l=cross(ol,s.l);
}
static St axpy(const St& a,const St& d,double t){ return {a.n1+d.n1*t,a.n2+d.n2*t,a.l+d.l*t}; }
int main(int argc,char** argv){
  const char* names[6]={"C","A","U1","U2","U3","U4"}; const int N=argc>1?std::atoi(argv[1]):300; basin=argc>5?std::atoi(argv[5]):0;
  std::vector<double> Bs{0.0,0.01,0.03,0.1,0.3,1.0,3.0}; if(argc>2) Bs={std::atof(argv[2])};
  for(double B:Bs){ zeta=(g/2)*muB*B/hbar;
    const double maxrate=4*PI*7.5e9+4*zeta; const double dt=std::min(2e-12,0.04/maxrate); const long steps=(long)std::ceil(5*tau/dt);
    const int v0=argc>3?std::atoi(argv[3]):0, v1=argc>4?std::atoi(argv[4]):4;
    for(variant=v0;variant<v1;++variant){ double res[3],se[3];
      for(int J=0;J<3;++J){ std::mt19937_64 rng(542+J); std::normal_distribution<double> nd; double sum=0,sum2=0;
        const double ct=((J+0.5)*(J+0.5)-2.25-1)/3, st=std::sqrt(std::max(0.0,1-ct*ct));
        for(int k=0;k<N;++k){ double q[4]; double qn=0; for(double& v:q){v=nd(rng); qn+=v*v;} qn=std::sqrt(qn); for(double& v:q) v/=qn;
          const double w=q[0],a=q[1],b=q[2],cq=q[3];
          auto R=[&](V v)->V{ V u{a,b,cq}; V t=cross(u,v)*2; return v+t*w+cross(u,t); };
          St s{R({st,0,ct}),R({st,0,ct}),R({0,0,1.5})};
          const V Jax=s.l+s.n1; const double xJ=4*zeta/std::abs(EJ[J]);
          const bool onJ=(basin==2)||(basin==3&&xJ<1);
          if(basin&&B>0){ const V sh=onJ?Jax*(1/nrm(Jax)):s.n1; const int m=(int)std::lround(1.5*sh[2]); V t;
            if(m==0){ t={sh[0],sh[1],0}; t=t*(1/nrm(t)); } else t={0,0,(double)m};
            const V ax=cross(sh,t); const double sa=nrm(ax), ca=dot(sh,t);
            if(sa>1e-12){ const V u=ax*(1/sa); auto Rot=[&](V v)->V{ return v*ca+cross(u,v)*sa+u*(dot(u,v)*(1-ca)); };
              s={Rot(s.n1),Rot(s.n2),Rot(s.l)}; } }
          double acc=0; St d1,d2,d3,d4;
          for(long i=0;i<steps;++i){ const double t=(i+0.5)*dt;
            deriv(s,d1); deriv(axpy(s,d1,dt/2),d2); deriv(axpy(s,d2,dt/2),d3); deriv(axpy(s,d3,dt),d4);
            St n=s; for(int c=0;c<3;++c){ n.n1[c]+=dt/6*(d1.n1[c]+2*d2.n1[c]+2*d3.n1[c]+d4.n1[c]); n.n2[c]+=dt/6*(d1.n2[c]+2*d2.n2[c]+2*d3.n2[c]+d4.n2[c]); n.l[c]+=dt/6*(d1.l[c]+2*d2.l[c]+2*d3.l[c]+d4.l[c]); }
            n.n1=n.n1*(1/nrm(n.n1)); n.n2=n.n2*(1/nrm(n.n2)); n.l=n.l*(1.5/nrm(n.l)); s=n;
            acc+=G*std::exp(-G*t)*dt*(1-dot(s.n1,s.n2))/2; }
          acc+=std::exp(-5.0)*(1-dot(s.n1,s.n2))/2;   // tail beyond 5 tau with the final w
          sum+=acc; sum2+=acc*acc; }
        res[J]=sum/N; se[J]=std::sqrt(std::max(0.0,sum2/N-res[J]*res[J])/N); }
      const double tot=(res[0]+3*res[1]+5*res[2])/9, tse=std::sqrt(se[0]*se[0]+9*se[1]*se[1]+25*se[2]*se[2])/9;
      std::printf("B %5.2f T %-2s%s: P_S J0 %.5e(%.1e) J1 %.5e(%.1e) J2 %.5e(%.1e) 3P %.5e(%.1e)\n",B,names[variant],basin==1?"+basin(s)":basin==2?"+basin(J)":basin==3?"+basin(x)":"",res[0],se[0],res[1],se[1],res[2],se[2],tot,tse);
      std::fflush(stdout); } }
}

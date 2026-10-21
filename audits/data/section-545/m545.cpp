// Audit 545 (from m543.cpp): Born projection on field eigenstates (eig545_<B>.txt) + adiabatic field ramp from the B = 0 label.
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
  else if((variant==4||variant==5)&&nrm(s.n1+s.n2)<2e-6){ }   // singlet branch (534) as in the engine for s <= 1e-6
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

static void step(St& s,double dt){ St d1,d2,d3,d4; deriv(s,d1); deriv(axpy(s,d1,dt/2),d2); deriv(axpy(s,d2,dt/2),d3); deriv(axpy(s,d3,dt),d4);
  St n=s; for(int c=0;c<3;++c){ n.n1[c]+=dt/6*(d1.n1[c]+2*d2.n1[c]+2*d3.n1[c]+d4.n1[c]); n.n2[c]+=dt/6*(d1.n2[c]+2*d2.n2[c]+2*d3.n2[c]+d4.n2[c]); n.l[c]+=dt/6*(d1.l[c]+2*d2.l[c]+2*d3.l[c]+d4.l[c]); }
  n.n1=n.n1*(1/nrm(n.n1)); n.n2=n.n2*(1/nrm(n.n2)); n.l=n.l*(1.5/nrm(n.l)); s=n; }
int main(int argc,char** argv){
  // argv: N B variant(4 = U3, 5 = U4) [Tramp ns]
  const int N=std::atoi(argv[1]); const double B=std::atof(argv[2]); variant=std::atoi(argv[3]); const double Tr=(argc>4?std::atof(argv[4]):20.0)*1e-9;
  char fn[64]; std::snprintf(fn,sizeof fn,"eig545_%g.txt",B); FILE* f=std::fopen(fn,"r"); if(!f){ std::printf("no %s\n",fn); return 1; }
  const double zetaB=(g/2)*muB*B/hbar; const double maxrate=4*PI*7.5e9+4*zetaB; const double dt=std::min(2e-12,0.04/maxrate);
  int kind,J,M; double t,E; double sumP=0,sumW=0,sumQ=0; int k=0;
  while(std::fscanf(f,"%d %d %d %lf %lf",&kind,&J,&M,&t,&E)==5){
    std::mt19937_64 rng(545+100*k); std::uniform_real_distribution<double> U(0,1); double acc=0,acc2=0;
    for(int r=0;r<N;++r){
      const double ph1=2*PI*U(rng), ph2=2*PI*U(rng), ph3=2*PI*U(rng);
      St s;
      if(kind==0){ // 3P_J M: J^ at cos a = M/J, S^ on its cone (cos(S,J) from the triangle), L = J - S
        const double Jl=J+0.5, ca=J>0?(double)M/J:2*U(rng)-1, sa=std::sqrt(std::max(0.0,1-ca*ca)); const V Jh{sa*std::cos(ph1),sa*std::sin(ph1),ca};   // basin convention: cos = M/J
        const double cb=(1+Jl*Jl-2.25)/(2*Jl), sb=std::sqrt(std::max(0.0,1-cb*cb));
        V e1=cross(Jh,V{0.3,0.5,0.81}); e1=e1*(1/nrm(e1)); const V e2=cross(Jh,e1);
        const V sh=Jh*cb+(e1*std::cos(ph2)+e2*std::sin(ph2))*sb; s={sh,sh,Jh*Jl-sh};
      } else {     // 1P_1 M: cos(L^, B) = M, antiparallel spins in a random direction
        const double cl=(double)M, sl=std::sqrt(std::max(0.0,1-cl*cl)); const V l{1.5*sl*std::cos(ph1),1.5*sl*std::sin(ph1),1.5*cl};
        const double cz=2*U(rng)-1, sz=std::sqrt(1-cz*cz); const V n{sz*std::cos(ph3),sz*std::sin(ph3),cz}; s={n,n*(-1.0),l};
      }
      const long nr=(long)std::ceil(Tr/dt); for(long i=0;i<nr;++i){ zeta=zetaB*(i+0.5)/nr; step(s,dt); }
      zeta=zetaB; const long steps=(long)std::ceil(5*tau/dt); double a=0;
      for(long i=0;i<steps;++i){ step(s,dt); a+=G*std::exp(-G*(i+1)*dt)*dt*(1-dot(s.n1,s.n2))/2; }
      a+=std::exp(-5.0)*(1-dot(s.n1,s.n2))/2; acc+=a; acc2+=a*a; }
    const double m=acc/N, se=std::sqrt(std::max(0.0,acc2/N-m*m)/N);
    std::printf("B %5.2f T %s k %2d %s%d M %+d: t %.4f s_QM %.4f w_cl %.4f(%.4f)\n",B,variant==4?"U3":"U4",k,kind?"1P":"3P",J,M,t,1-t,m,se);
    sumP+=t/9; sumW+=t/9*m; sumQ+=t/9*(1-t); ++k; }
  std::printf("B %5.2f T %s: 3P Born sum_k p_k w_cl %.5e, QM %.5e, ratio %.3f (sum p %.4f, T_ramp %.0f ns)\n",B,variant==4?"U3":"U4",sumW,sumQ,sumW/sumQ,sumP,Tr*1e9);
}
